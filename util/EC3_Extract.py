import argparse

import openpyxl as xl

import numpy as np

parser = argparse.ArgumentParser(description="Extract section parameters from Blue Book spreadsheet.")

parser.add_argument('--section',      help='Specify section type.',    default=None, choices=['UB', 'RHS', 'SHS', 'CHS'])
parser.add_argument('--source',       help='Specify input file name.', type=str)
parser.add_argument('--name',         help='Specify internal name.',   type=str, default='unspecified')

args = parser.parse_args()

wb = xl.load_workbook(filename=args.source)
ws = wb.active

dataBegin = False
dataCount = 0
dataMajor = None

def add_minor_UB(row, major):
    B = float(row[ 5]) / 1E3
    D = float(row[ 4]) / 1E3
    F = float(row[ 7]) / 1E3
    W = float(row[ 6]) / 1E3
    R = float(row[ 8]) / 1E3
    Y = float(row[17]) / 1E8
    Z = float(row[18]) / 1E8
    J = float(row[28]) / 1E8
    A = float(row[29]) / 1E4
    print('        S[\'{major}\'][\'{minor}\'] = Universal({b:.3e}, {d:.3e}, {f:.3e}, {w:.3e}, {r:.3e}, area={a:.3e}, Iyy={y:.3e}, Izz={z:.3e}, J={j:.3e})'.format(major=major, minor=row[1], b=B, d=D, f=F, w=W, r=R, a=A, y=Y, z=Z, j=J))

def add_minor_RHS(row, major):
    wh = major.split()
    W = float(wh[0]) / 1E3
    H = float(wh[2]) / 1E3
    T = float(row[ 1]) / 1E3
    A = float(row[ 4]) / 1E4
    Y = float(row[ 7]) / 1E8
    Z = float(row[ 8]) / 1E8
    J = float(row[15]) / 1E8
    R = 2 * T
    print('        S[\'{major}\'][\'{minor}\'] = RHS({w:.3e}, {h:.3e}, {t:.3e}, {r:.3e}, area={a:.3e}, Iyy={y:.3e}, Izz={z:.3e}, J={j:.3e})'.format(major=major, minor=row[1], w=W, h=H, t=T, r=R, a=A, y=Y, z=Z, j=J))

def add_minor_SHS(row, major):
    wh = major.split()
    W = float(wh[0]) / 1E3
    H = float(wh[2]) / 1E3
    T = float(row[ 1]) / 1E3
    A = float(row[ 4]) / 1E4
    I = float(row[ 6]) / 1E8
    J = float(row[10]) / 1E8
    R = 2 * T
    print('        S[\'{major}\'][\'{minor}\'] = RHS({w:.3e}, {h:.3e}, {t:.3e}, {r:.3e}, area={a:.3e}, Iyy={y:.3e}, Izz={z:.3e}, J={j:.3e})'.format(major=major, minor=row[1], w=W, h=H, t=T, r=R, a=A, y=I, z=I, j=J))

def add_minor_CHS(row, major):
    D = float(major) / 1E3
    T = float(row[ 1]) / 1E3
    A = float(row[ 4]) / 1E4
    I = float(row[ 6]) / 1E8
    J = float(row[10]) / 1E8
    print('        S[\'{major}\'][\'{minor}\'] = CHS({d:.3e}, {t:.3e}, area={a:.3e}, Iyy={y:.3e}, Izz={z:.3e}, J={j:.3e})'.format(major=major, minor=row[1], d=D, t=T, a=A, y=I, z=I, j=J))

for row in ws.values:
    if not dataBegin:
        if args.section == 'UB':
            if row[4] == 'mm':
                dataBegin = True
        elif row[0] == 'mm':
            dataBegin = True

        if dataBegin:
            print('    __' + args.name + ' = None')
            print('')
            print('    @staticmethod')
            print('    def __create_' + args.name + '():')
            print('        S = {}')
            continue # starts with next row

    if not dataBegin:
        continue # still searching

    if row[0] is None: # no more data
        print('')
        print('        return S # No. data rows = {c}'.format(c=dataCount))
        print('')
        print('    @staticmethod')
        print('    def ' + args.name + '():')
        print('        if BlueBook.__' + args.name + ' is None:')
        print('            BlueBook.__' + args.name + ' = BlueBook.__create_' + args.name + '()')
        print('        return BlueBook.__' + args.name)
        print('')
        break

    dataCount += 1

    if row[0] != '':
        dataMajor = row[0]
        print('        S[\'' + dataMajor + '\'] = {}')

    if args.section == 'UB':
        add_minor_UB(row, dataMajor)
    elif args.section == 'RHS':
        add_minor_RHS(row, dataMajor)
    elif args.section == 'SHS':
        add_minor_SHS(row, dataMajor)
    elif args.section == 'CHS':
        add_minor_CHS(row, dataMajor)
