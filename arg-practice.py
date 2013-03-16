#!/usr/bin/env python

import argparse

parser = argparse.ArgumentParser()

parser.add_argument( 't_ty', help="hello helpo" )
parser.add_argument( '-n', '--numm', type=int, help="numero" )

m_group = parser.add_mutually_exclusive_group()
m_group.add_argument( '--poo', default="who?" )
a_group = m_group.add_argument_group()
a_group.add_argument( '--kelly', help="kelly!" )
a_group.add_argument( '--res1', type=int )
a_group.add_argument( '--res2', type=int )


args = parser.parse_args()

print args.poo
print args.t_ty
if args.numm:
    print args.numm + 2

if args.kelly:
    print args.res1
    print args.res2
