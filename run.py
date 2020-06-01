from argparse import ArgumentParser
import bch


def run():
    '''Parse argument of the command line'''
    parser = ArgumentParser(prog='coder', description='Bose–Chaudhuri–Hocquenghem codes')
    subparsers = parser.add_subparsers(title='mode', description='Working mode')

    generation_mode = subparsers.add_parser('generate',
                                            help='Generates BCH code with given parameters')
    generation_mode.add_argument('-n', dest='n', type=int, help='desired maximum length of a '
                                 + 'message block transmitted through a communication channel',
                                 required=True)
    generation_mode.add_argument('-p', dest='p', type=float, help='probability of an error in the '
                                 + 'communication channel', required=True)
    generation_mode.add_argument('-o', dest='out_file', type=str, default='code.data',
                                 help='file with information for coder and decoder', required=True)
    generation_mode.set_defaults(func=bch.generate)

    encode_mode = subparsers.add_parser('encode', help='Encodes message with an error')
    encode_mode.add_argument('-i', dest='in_file', type=str, help='file with information for coder',
                             required=True)
    encode_mode.add_argument('-o', dest='out_file', type=str, default='code.data',
                             help='output file consists encoded message with error', required=True)
    encode_mode.set_defaults(func=bch.encode)

    decode_mode = subparsers.add_parser('decode', help='Finds error and decodes message')
    decode_mode.add_argument('-i', dest='in_file', type=str, help='file with information for decoder',
                             required=True)
    decode_mode.add_argument('-m', dest='message_file', type=str, help='file consists message to decode',
                             required=True)
    decode_mode.set_defaults(func=bch.decode)

    args = parser.parse_args()
    args.func(args)

if __name__ == "__main__":
    run()
