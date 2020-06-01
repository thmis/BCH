import numpy as np
import math
from finite_field import GaloisField
import math_utils
import pickle
import random

from itertools import combinations

class BCH:
    '''Generates gererator polynomial for a BCH code
       with given parameters

       Parameters:
        n — desired maximum length of a message block transmitted
            through a communication channel
        p — probability of an error in the communication channel'''

    def initialize(self, n: int, t: int, p: float, power: int, k: int, generator):
        self.n = n
        self.t = t
        self.gf = GaloisField(power)
        self.k = k
        self.generator = generator
        self.p = p


    def generate(self, n: int, p: float):
        power = math.floor(math.log2(n + 1))
        self.p = p
        self.n = 2 ** power - 1
        while self.n * p > power - 1:
            power -= 1
            self.n = 2 ** power - 1

        self.power = power
        self.t = math.ceil(self.n * self.p)
        self.k = self.n - self.t * power
        self.gf = GaloisField(power)
        self.generator = self.calculate_generator_polynomial()


    def write_to_file(self, file_name: str):
        '''Write information about the generated code to the file'''
        code_json = {}
        code_json['n'] = self.n
        code_json['t'] = self.t
        code_json['k'] = self.k
        code_json['p'] = self.p
        code_json['power'] = self.power
        code_json['Generator'] = self.generator

        with open(file_name, 'wb') as file:
            pickle.dump(code_json, file)


    def calculate_generator_polynomial(self):
        generator = self.gf.primitive_polynomial
        size = math_utils.greatest_bit(generator) - 1
        for i in range(2, self.t + 1):
            min_pol = self.gf.get_minimal_polynomial(self.gf.cyclotomic_classes[i])
            generator = math_utils.multiply_polynomials(generator, min_pol)

        return generator[-math_utils.greatest_bit(generator):]


    def str_to_bits(self, s):
        result = []
        for c in s:
            bits = bin(ord(c))[2:]
            bits = '00000000'[len(bits):] + bits
            result.extend([int(b) for b in bits])
        return result


    def bits_to_str(self, bits):
        chars = []
        for b in range(int(len(bits) / 8)):
            byte = bits[b * 8 : (b+1) * 8]
            chars.append(chr(int(''.join([str(bit) for bit in byte]), 2)))
        return ''.join(chars)


    def find_cycl_class(self, num):
        for i in range(len(self.gf.cyclotomic_classes)):
            if num in self.gf.cyclotomic_classes[i]:
                return i
        return -1


    def get_syndromes(self, message):
        syndromes = []
        for i in range(1, 2 * self.t + 1):
            min_pol = self.gf.get_minimal_polynomial(self.gf.cyclotomic_classes[self.find_cycl_class(i)])
            syndrome_pol = math_utils.divide_polynomials(message, min_pol)[1]
            syndrome = math_utils.get_powered_polynomial(syndrome_pol, i)
            syndrome = math_utils.divide_polynomials(syndrome, self.gf.primitive_polynomial)[1]
            syndromes.append(math_utils.bin_vector_to_num(syndrome))
        syndromes = [self.gf.reversed_log_table[syndrome] for syndrome in syndromes]
        return syndromes


    def get_locators_polynomial(self, syndromes):
        sigma = np.zeros(2 * self.t - 1, dtype=int)
        sigma = np.concatenate((np.array([1]), sigma))
        counter = 0
        B = sigma.copy()

        for j in range(2 * self.t):
            b = 0
            for i in range(0, counter + 1):
                if not (sigma[i] and math_utils.bin_vector_to_num(
                    self.gf.logarithm_table[syndromes[j - i]])):
                    continue
                
                sigma_i = self.gf.reversed_log_table[sigma[i]]
                mult = (sigma_i + syndromes[j - i]) % self.gf.field_power
                b ^= math_utils.bin_vector_to_num(self.gf.logarithm_table[mult])

            for i in range(len(B) - 1, -1, -1):
                B[i] = 0 if (i == 0) else B[i - 1]
            if b != 0:
                temp = sigma.copy()
                for i in range(len(sigma)):
                    B_i = self.gf.reversed_log_table[B[i]]
                    b_idx = self.gf.reversed_log_table[b]
                    if B_i != -1:
                        mult = (b_idx + B_i) % self.gf.field_power
                        temp[i] ^= math_utils.bin_vector_to_num(self.gf.logarithm_table[mult])
                if 2*counter <= j:
                    for i in range(len(sigma)):
                        B[i] = 0
                        delta_idx = self.gf.field_power - self.gf.reversed_log_table[b]
                        sigma_i = self.gf.reversed_log_table[sigma[i]]
                        if sigma_i != -1:
                            div = (sigma_i + delta_idx) % self.gf.field_power
                            B[i] = math_utils.bin_vector_to_num(self.gf.logarithm_table[div])
                    sigma = temp.copy()
                    counter = j + 1 - counter
                else:
                    sigma = temp.copy()
        result = []
        for poly in sigma:
            result.append(self.gf.reversed_log_table[poly])
        return result


    def decode_block(self, block):
        decode_succeeded = True
        syndromes = self.get_syndromes(block)
        decoded = block
        locator_polynomial = self.get_locators_polynomial(syndromes)
        for i in range(len(locator_polynomial)):
            if locator_polynomial[i] == -1:
                locator_polynomial = locator_polynomial[:i]
                break
        roots = self.gf.find_roots_of_polynomial(locator_polynomial)
        if sum(locator_polynomial) != 0 and len(roots) == 0:
            decode_succeeded = False
        locations = []
        for root in roots:
            locations.append((self.gf.field_power - root) % self.gf.field_power)
        for location in locations:
            decoded[-location-1] ^= 1
        return decoded[:self.k], decode_succeeded


    def encode_block(self, block):
        power = math_utils.greatest_bit(self.generator)
        encoded = np.concatenate((block, np.zeros(power - 1, dtype=int)))
        encoded = math_utils.xor(encoded, math_utils.divide_polynomials(encoded, self.generator)[1],
                len(encoded))
        return encoded

        
def construct_code_from_file(file_name: str):
    with open(file_name, 'rb') as file:
        code_json = pickle.load(file)

    n = code_json['n']
    t = code_json['t']
    k = code_json['k']
    p = code_json['p']
    power = code_json['power']
    generator = code_json['Generator']

    code = BCH()
    code.initialize(n, t, p, power, k, generator)
    return code


def encode(args):
    random.seed()
    code = construct_code_from_file(args.in_file)
    message = input('Input the message for encoding: \n >')
    encoded_blocks = []

    binary = code.str_to_bits(message)
    while len(binary) % code.k != 0:
        binary.append(0)

    for i in range(int(len(binary) / code.k)):
        src_line = binary[i * code.k: (i+1) * code.k]
        encoded = code.encode_block(src_line)
        encoded_blocks.append(encoded)

        print(src_line, list(encoded), sep=' ')

    with open(args.out_file, 'w') as out_file:
        for block in encoded_blocks:
            error_vector = np.zeros(code.n, dtype=int)
            for i in range(len(error_vector)):
                random_number = random.randint(1, 100)
                if random_number < 100*code.p:
                    error_vector[i] = 1
            block_with_errors = math_utils.xor(block, error_vector, code.n)
            print(block, error_vector, block_with_errors, sep=' ')
            for bit in block_with_errors:
                out_file.write(str(bit))
    

def decode(args):
    code = construct_code_from_file(args.in_file)
    with open(args.message_file, 'r') as in_file:
        message = in_file.read()
        decoded = np.array([], dtype=int)
        for i in range(int(len(message) / code.n)):
            src_line = message[i * code.n: (i+1) * code.n]
            block = np.array([], dtype=int)
            for ch in src_line:
                block = np.append(block, int(ch))
            
            decoded_block, success = code.decode_block(block)
            print(block, decoded_block, success, sep=' ')
            decoded = np.concatenate((decoded, decoded_block))


    plaintext = code.bits_to_str(decoded)
    print(plaintext)


def generate(args):
    code = BCH()
    code.generate(args.n, args.p)
    code.write_to_file(args.out_file)
