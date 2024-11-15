import random
import argparse


def random_DNA_gen(x):
    return ''.join(random.choice('GCAT') for _ in range(x))

def valid_barcode(barcode):
    return 3 <= sum(c in 'GC' for c in barcode) <= 5

def generate_unique_barcodes(num_barcodes, barcode_length):
    barcodes = set()
    while len(barcodes) < num_barcodes:
        barcode = random_DNA_gen(barcode_length)
        if barcode not in barcodes and valid_barcode(barcode):
            barcodes.add(barcode)
    return list(barcodes)

def main():
    parser = argparse.ArgumentParser(description="barcode generator")
    parser.add_argument("-bl", "--barcode_length", type=int, default=8,
                        help="Length of the barcode")
    parser.add_argument("-nb", "--num_barcodes", type=int, default=200,
                        help="Number of barcodes generated")
    args = parser.parse_args()

    unique_barcodes = generate_unique_barcodes(args.num_barcodes, args.barcode_length)
    return unique_barcodes

if __name__ == "__main__":
    barcodes = main()
    with open('barcodes.txt', 'w') as file:
        for barcode in barcodes:
            file.write(barcode + '\n')
