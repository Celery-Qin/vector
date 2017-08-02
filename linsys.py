# coding=utf-8

import pdb
from decimal import Decimal, getcontext
from vector import Vector
from plane import Plane
from copy import deepcopy

getcontext().prec = 30


class LinearSystem(object):

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = \
        'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:  # 确保每一个平面的维度都和第一个平面的维度相等
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)

    def swap_rows(self, row1, row2):
        """将row1和row2位置调换"""
        self[row1], self[row2] = self[row2], self[row1]

    def multiply_coefficient_and_row(self, coefficient, row):
        """将row乘以coefficient，新建一个plane，若在老的plane上修改，还要重新初始化basepoint"""
        n = self[row].normal_vector
        k = self[row].constant_term

        new_normal_vector = n.scalar_multiply(coefficient)
        new_constant_term = k * coefficient
        self[row] = Plane(new_normal_vector, new_constant_term)

    def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):
        """将row_to_add乘以coefficient后，加到row_to_be_added_to上"""
        n1 = self[row_to_add].normal_vector
        n2 = self[row_to_be_added_to].normal_vector
        k1 = self[row_to_add].constant_term
        k2 = self[row_to_be_added_to].constant_term

        new_normal_vector = n1.scalar_multiply(coefficient).plus(n2)
        new_constant_term = (k1 * coefficient) + k2

        self[row_to_be_added_to] = Plane(new_normal_vector, new_constant_term)

    def indices_of_first_nonzero_terms_in_each_row(self):
        """求每行首个不为零系数项的index"""
        num_equations = len(self)  # 有几个平面
        num_variables = self.dimension  # 有几个变量（维度）

        indices = [-1] * num_equations  # 设为[-1,-1,-1]

        for i, p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices

    def __len__(self):
        return len(self.planes)

    def __getitem__(self, i):
        return self.planes[i]

    def __setitem__(self, i, x):
        """把self中第i个平面换成x"""
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)

    def __str__(self):
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i + 1, p)
                for i, p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret

    def compute_triangular_form(self):
        """
           1、只能与最接近的row交换顺序
           2、不要只将row与数相乘
           3、只将几倍的row与下面的row相加
           A*x_1 + B*x_2 + C*x_3 = k_1 
           D*x_1 + E*x_2 + F*x_3 = k_2
           G*x_1 + H*x_2 + I*x_3 = k_3
        """
        system = deepcopy(self)  # 用copy，不对原始方程组做修改，以防以后还要用到它

        num_equations = len(system)  # 有几条方程
        num_variables = system.dimension  # 有几个变量

        j = 0
        for i in range(num_equations):
            while j < num_variables:
                c = MyDecimal(system[i].normal_vector[j])  # 第i条方程，第j个变量
                if c.is_near_zero():
                    swap_succeeded = system.swap_with_row_below_for_nonzero_coefficient_if_able(
                        i, j)
                    if not swap_succeeded:  # 如果没有交换成功，说明下面方程的这个variable系数也都是0
                        j += 1  # 就不看这个系数啦，跳过
                        continue

                system.clear_coefficient_below(i, j)
                j += 1  # 继续看下一个系数
                break

        return system

    def swap_with_row_below_for_nonzero_coefficient_if_able(self, row, col):
        num_equations = len(self)

        for k in range(row + 1, num_equations):
            coefficient = MyDecimal(self[k].normal_vector[col])
            if not coefficient.is_near_zero():
                self.swap_rows(row, k)
                return True
        return False

    def clear_coefficient_below(self, row, col):
        num_equations = len(self)
        beta = MyDecimal(self[row].normal_vector[col])

        for k in range(row + 1, num_equations):
            n = self[k].normal_vector
            gamma = n[col]
            alpha = -gamma / beta  # 当前行乘以alpha，并与下面的行相加
            self.add_multiple_times_row_to_row(alpha, row, k)

    def compute_rref(self):
        """
        Reduce Row-Echelon Form:
        -Triangular form
        -Each pivot variable has coefficient 1
        -Each pivot variable is in own column
        -(Technical)Any 0=k equations should be 0=1
		 should come before 0=0
        """
        return deepcopy(self)


class MyDecimal(Decimal):

    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps

# 测试row运算函数
# p0 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='1')
# p1 = Plane(normal_vector=Vector(['0', '1', '0']), constant_term='2')
# p2 = Plane(normal_vector=Vector(['1', '1', '-1']), constant_term='3')
# p3 = Plane(normal_vector=Vector(['1', '0', '-2']), constant_term='2')

# s = LinearSystem([p0, p1, p2, p3])

# s.swap_rows(0, 1)
# # pdb.set_trace()
# # print s[1] == p0
# if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
#     print 'test case 1 failed'

# s.swap_rows(1, 3)
# if not (s[0] == p1 and s[1] == p3 and s[2] == p2 and s[3] == p0):
#     print 'test case 2 failed'

# s.swap_rows(3, 1)
# if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
#     print 'test case 3 failed'

# s.multiply_coefficient_and_row(1, 0)
# if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
#     print 'test case 4 failed'

# s.multiply_coefficient_and_row(-1, 2)
# if not (s[0] == p1 and
#         s[1] == p0 and
#         s[2] == Plane(normal_vector=Vector(['-1', '-1', '1']), constant_term='-3') and
#         s[3] == p3):
#     print 'test case 5 failed'

# s.multiply_coefficient_and_row(10, 1)
# if not (s[0] == p1 and
#         s[1] == Plane(normal_vector=Vector(['10', '10', '10']), constant_term='10') and
#         s[2] == Plane(normal_vector=Vector(['-1', '-1', '1']), constant_term='-3') and
#         s[3] == p3):
#     print 'test case 6 failed'

# s.add_multiple_times_row_to_row(0, 0, 1)
# if not (s[0] == p1 and
#         s[1] == Plane(normal_vector=Vector(['10', '10', '10']), constant_term='10') and
#         s[2] == Plane(normal_vector=Vector(['-1', '-1', '1']), constant_term='-3') and
#         s[3] == p3):
#     print 'test case 7 failed'

# s.add_multiple_times_row_to_row(1, 0, 1)
# if not (s[0] == p1 and
#         s[1] == Plane(normal_vector=Vector(['10', '11', '10']), constant_term='12') and
#         s[2] == Plane(normal_vector=Vector(['-1', '-1', '1']), constant_term='-3') and
#         s[3] == p3):
#     print 'test case 8 failed'

# s.add_multiple_times_row_to_row(-1, 1, 0)
# if not (s[0] == Plane(normal_vector=Vector(['-10', '-10', '-10']), constant_term='-10') and
#         s[1] == Plane(normal_vector=Vector(['10', '11', '10']), constant_term='12') and
#         s[2] == Plane(normal_vector=Vector(['-1', '-1', '1']), constant_term='-3') and
#         s[3] == p3):
#     print 'test case 9 failed'

# 三角形函数--变形！
# p1 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='1')
# p2 = Plane(normal_vector=Vector(['0', '1', '1']), constant_term='2')
# s = LinearSystem([p1, p2])
# t = s.compute_triangular_form()

# if not (t[0] == p1 and
#         t[1] == p2):
#     print 'test case 1 failed'

# p1 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='1')
# p2 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='2')
# s = LinearSystem([p1, p2])
# t = s.compute_triangular_form()
# # pdb.set_trace()
# if not (t[0] == p1 and
#         t[1] == Plane(constant_term='1')):
#     print 'test case 2 failed'

# p1 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='1')
# p2 = Plane(normal_vector=Vector(['0', '1', '0']), constant_term='2')
# p3 = Plane(normal_vector=Vector(['1', '1', '-1']), constant_term='3')
# p4 = Plane(normal_vector=Vector(['1', '0', '-2']), constant_term='2')
# s = LinearSystem([p1, p2, p3, p4])
# t = s.compute_triangular_form()
# if not (t[0] == p1 and
#         t[1] == p2 and
#         t[2] == Plane(normal_vector=Vector(['0', '0', '-2']), constant_term='2') and
#         t[3] == Plane()):
#     print 'test case 3 failed'

# p1 = Plane(normal_vector=Vector(['0', '1', '1']), constant_term='1')
# p2 = Plane(normal_vector=Vector(['1', '-1', '1']), constant_term='2')
# p3 = Plane(normal_vector=Vector(['1', '2', '-5']), constant_term='3')
# s = LinearSystem([p1, p2, p3])
# t = s.compute_triangular_form()
# if not (t[0] == Plane(normal_vector=Vector(['1', '-1', '1']), constant_term='2') and
#         t[1] == Plane(normal_vector=Vector(['0', '1', '1']), constant_term='1') and
#         t[2] == Plane(normal_vector=Vector(['0', '0', '-9']), constant_term='-2')):
#     print 'test case 4 failed'

# 编写REFF函数
p1 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['0', '1', '1']), constant_term='2')
s = LinearSystem([p1, p2])
r = s.compute_rref()
# pdb.set_trace()
if not (r[0] == Plane(normal_vector=Vector(['1', '0', '0']), constant_term='-1') and
        r[1] == p2):
    print 'test case 1 failed'

p1 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='2')
s = LinearSystem([p1, p2])
r = s.compute_rref()
if not (r[0] == p1 and
        r[1] == Plane(constant_term='1')):
    print 'test case 2 failed'

p1 = Plane(normal_vector=Vector(['1', '1', '1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['0', '1', '0']), constant_term='2')
p3 = Plane(normal_vector=Vector(['1', '1', '-1']), constant_term='3')
p4 = Plane(normal_vector=Vector(['1', '0', '-2']), constant_term='2')
s = LinearSystem([p1, p2, p3, p4])
r = s.compute_rref()
if not (r[0] == Plane(normal_vector=Vector(['1', '0', '0']), constant_term='0') and
        r[1] == p2 and
        r[2] == Plane(normal_vector=Vector(['0', '0', '-2']), constant_term='2') and
        r[3] == Plane()):
    print 'test case 3 failed'

p1 = Plane(normal_vector=Vector(['0', '1', '1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['1', '-1', '1']), constant_term='2')
p3 = Plane(normal_vector=Vector(['1', '2', '-5']), constant_term='3')
s = LinearSystem([p1, p2, p3])
r = s.compute_rref()
if not (r[0] == Plane(normal_vector=Vector(['1', '0', '0']), constant_term=Decimal('23') / Decimal('9')) and
        r[1] == Plane(normal_vector=Vector(['0', '1', '0']), constant_term=Decimal('7') / Decimal('9')) and
        r[2] == Plane(normal_vector=Vector(['0', '0', '1']), constant_term=Decimal('2') / Decimal('9'))):
    print 'test case 4 failed'
