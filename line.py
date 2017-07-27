# coding=utf-8

from decimal import Decimal, getcontext

from vector import Vector

getcontext().prec = 30


class Line(object):
    """docstring for Line"""

    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'  # 找不到非零元素（全是零）的警告

    def __init__(self, normal_vector=None, constant_term=None):
        """输入直线法向量、直线等式常量（标准式等号右侧k）。
        法向量给出了标准直线形式的系数。"""
        self.dimension = 2  # 只考虑二维直线

        if not normal_vector:  # 如果法向量为空，设标准式系数皆为0，那么这个直线的法向量是零向量
            all_zeros = ['0'] * self.dimension
            normal_vector = Vector(all_zeros)
        self.normal_vector = normal_vector

        if not constant_term:  # 如果标准式常量为空，设为0，并设置精确度
            constant_term = Decimal('0')
        self.constant_term = Decimal(constant_term)

        self.set_basepoint()

    def is_parallel_to(self, ell):
        n1 = self.normal_vector
        n2 = ell.normal_vector
        return n1.parallelism(n2)

    def set_basepoint(self):
        """设置直线的基准点，方向向量和基准点决定一条直线"""
        try:
            n = self.normal_vector
            c = self.constant_term
            basepoint_coords = ['0'] * self.dimension  # 默认基准点为原点

            initial_index = Line.first_nonzero_index(n)  # 找到标准式中第一个不为0系数的序号
            initial_coefficient = n[initial_index]  # 标准式中第一个不为0的系数

            basepoint_coords[initial_index] = c / initial_coefficient
            self.basepoint = Vector(basepoint_coords)

        except Exception as e:
            if str(e) == Line.NO_NONZERO_ELTS_FOUND_MSG:
                self.basepoint = None
            else:
                raise e

    def __str__(self):
        """写出直线标准式"""
        num_decimal_places = 3

        def write_coefficient(coefficient, is_initial_term=False):
            coefficient = round(coefficient, num_decimal_places)  # 系数四舍五入到三位
            if coefficient % 1 == 0:
                coefficient = int(coefficient)

            output = ''

            if coefficient < 0:
                output += '-'
            if coefficient > 0 and not is_initial_term:
                output += '+'

            if not is_initial_term:
                output += ' '
            if abs(coefficient) != 1:  # 如果系数的绝对值为1，是不用写出来的
                output += '{}'.format(abs(coefficient))

            return output

        n = self.normal_vector

        try:
            initial_index = Line.first_nonzero_index(n)  # 第一个不为零的系数的索引
            terms = [write_coefficient(n[i], is_initial_term=(i == initial_index)) + 'x_{}'.format(i + 1)
                     for i in range(self.dimension) if round(n[i], num_decimal_places) != 0]
            output = ' '.join(terms)

        except Exception as e:
            if str(e) == self.NO_NONZERO_ELTS_FOUND_MSG:
                output = '0'
            else:
                raise e
        constant = round(self.constant_term, num_decimal_places)
        if constant % 1 == 0:
            constant = int(constant)
        output += '= {}'.format(constant)

        return output

    def __eq__(self, ell):
        if self.is_parallel_to(ell) and self.basepoint.__eq__(ell.basepoint):
            return True
        else:
            return False

    def intersection_with(self, ell):
        pass

    @staticmethod
    def first_nonzero_index(iterable):
        """找到等式的第一个非零系数"""
        for k, item in enumerate(iterable):
            if not MyDecimal(item).is_near_zero():  # 如果item不为零，则返回它的序号
                return k
        raise Exception(Line.NO_NONZERO_ELTS_FOUND_MSG)


class MyDecimal(Decimal):
    """检测某个数字是否在误差范围0内，避免因为四舍五入错误而出现错误的答案"""

    def is_near_zero(self, eps=1e-10):  # 如果一个值的绝对值小于1e-10，认为它near zero
        return abs(self) < eps

# Q1:Determine if two lines are parallel
# Q2:Determine if two lines are equal
# Q3:Compute the intersection of two lines
# Or return some indication of no intersection/infinite intersection

Line1 = Line(Vector([4.046, 2.836]), 1.21)
Line2 = Line(Vector([10.115, 7.09]), 3.025)
print 'Are they parallel?', Line1.is_parallel_to(Line2)
print 'Are they equal?', Line1.__eq__(Line2)
print 'the intersection is ', Line1.intersection_with(Line2)
