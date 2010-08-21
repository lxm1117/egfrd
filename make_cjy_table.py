#!/usr/bin/python


import math
import numpy
import scipy.special as special

import string



def minz_j(n):
    
    # there is no singularity in besselJ at z=0
    return 0

def minz_y(n):
    # from GSL 1.12 bessel_Yn.c:
    # if(x < 5.0) {
    #   int status = bessel_Yn_small_x(n, x, result);
    return 5.

def maxz_j(n):

    # from GSL 1.12 bessel_Jn.c
    # else if(GSL_ROOT4_DBL_EPSILON * x > (n*n+1.0)) {
    #  int status = gsl_sf_bessel_Jnu_asympx_e((double)n, x, result);
    #  ...

    z = (n * n + 1) / 1.221e-4

    # but z can be too large..
    if z >= 1000:
        z = max(1000, n * n)

    return z


def maxz_y(n):

    # from GSL 1.12 bessel_Yn.c
    #  else if(GSL_ROOT3_DBL_EPSILON * x > (l*l + l + 1.0)) {
    #     int status = gsl_sf_bessel_Ynu_asympx_e(l + 0.5, x, result);
      #     ...
    z = (n * n + n + 1) / 6.06e-6

    # ... but this is usually too big.
    if z >= 2000:
        z = max(2000, n * n)

    return z


def JnYn(maxn, resolution):

    delta = numpy.pi / resolution
    z_table = numpy.mgrid[min(minz_j(maxn),minz_y(maxn)):
                              max(maxz_j(maxn), maxz_y(maxn)):delta]

    J_table = []
    Jdot_table = []
    Y_table = []
    Ydot_table = []
    for n in range(maxn+1):
        J_table.append(special.jn(n, z_table)) 
        Jdot_table.append(special.jvp(n, z_table))
        Y_table.append(special.yn(n, z_table))
        Ydot_table.append(special.yvp(n, z_table))

    return z_table, J_table, Jdot_table, Y_table, Ydot_table

def write_header(file):

    template = '''#ifndef CYLINDRICAL_BESSEL_TABLE_HPP
#define CYLINDRICAL_BESSEL_TABLE_HPP

/* Auto-generated by a script.  Do not edit. */

namespace cb_table
{

struct Table
{
    const unsigned int N;
    const double x_start;
    const double delta_x;
    const double* const y;
};
'''

    file.write(template)

def write_footer(file):

    template = '''
}  // namespace cbjy_table

#endif /* CYLINDRICAL_BESSEL_TABLE_HPP */
'''

    file.write(template)

def write_table_array(file, name, minn, maxn):

    file.write('static unsigned int %s_min(%d);\n' % (name, minn)) 
    file.write('static unsigned int %s_max(%d);\n' % (name, maxn)) 
    file.write('static const Table* %s[%d + 1] =\n{\n' % (name, maxn)) 

    for n in range(minn):
        file.write('    0,\n')

    for n in range(minn, maxn+1):
        file.write('    &%s%d,\n' % (name, n))

    file.write('};\n\n')


def write_array(file, name, table):

    head_template = '''
static const double %s[%d + 1] =
{\n'''

    number_template = '''    %.18f'''
    foot_template = '''};\n'''

    N = len(table)

    file.write(head_template % (name, N))
    file.write(',\n'.join([number_template % n for n in table]))
    file.write(foot_template)

def write_arrays(file, name, table1, table2):

    head_template = '''
static const double %s[%d + 1] =
{\n'''

    number_template = '''    %.18e, %.18e'''
    foot_template = '''};\n'''

    # check if len(table1) == len(table2)
    N = len(table1)

    file.write(head_template % (name, N * 2))

    file.write(',\n'.join([number_template % (value, table2[i]) for i, value in enumerate(table1)]))

    file.write(foot_template)


def write_table(file, name, N, x_start, delta_x):

    struct_template = '''
static const Table %s = { %d, %.18f, %.18f, %s_f };

'''

    file.write(struct_template % (name, N, x_start, delta_x, name))



if __name__ == '__main__':

    import sys

    filename = sys.argv[1]

    file = open(filename, 'w')

    minn_j = 0
    # this should be larger (than maxn_y), but the table bloats.
    maxn_j = 50

    minn_y = 0

    maxn_y = 50


    resolution = 35

    write_header(file)

    z_table, j_table, jdot_table, y_table, ydot_table \
        = JnYn(max(maxn_j, maxn_y), resolution)

    delta_z = z_table[1]-z_table[0]

    # j
    for n in range(minn_j, maxn_j + 1):

        start = numpy.searchsorted(z_table, minz_j(n))
        end = numpy.searchsorted(z_table, maxz_j(n))
        z_start = z_table[start]
        j = j_table[n][start:end]
        jdot = jdot_table[n][start:end]
        write_arrays(file, 'cj_table%d_f' % n, j, jdot)
        write_table(file, 'cj_table%d' % n, end-start, z_start, delta_z)
        print 'j', n, len(j)

    # y
    for n in range(minn_y, maxn_y + 1):

        start = numpy.searchsorted(z_table, minz_y(n))
        end = numpy.searchsorted(z_table, maxz_y(n))
        z_start = z_table[start]
        y = y_table[n][start:end]
        ydot = ydot_table[n][start:end]
        write_arrays(file, 'cy_table%d_f' % n, y, ydot)
        write_table(file, 'cy_table%d' % n, end-start, z_start, delta_z)

        print 'y', n, len(y)

    write_table_array(file, 'cj_table', minn_j, maxn_j)
    write_table_array(file, 'cy_table', minn_y, maxn_y)

    write_footer(file)

    file.write('\n')

    file.close()

