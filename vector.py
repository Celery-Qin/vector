# coding=utf-8

import pdb
import math
from decimal import Decimal, getcontext

getcontext().prec = 15


class Vector(object):

    CANNOT_NORMALIZE_ZERO_VECTOR_MSG = 'Cannot normalize the zero vector'
    NO_UNIQUE_PARALLEL_COMPONENT_MEG = 'No unique parallel component to the zero vector'
    ONLY_DEFINE_IN_TWO_THREE_DIMS_MSG = 'only defined in 2 and 3 dimensions'

    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = list([Decimal(x) for x in coordinates])
            self.dimension = len(coordinates)
            # self.idx = 0

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')

    def __iter__(self):
        self.current = -1
        return self

    def next(self):
        self.current += 1
        # self.idx += 1
        if self.current >= self.dimension:
            raise StopIteration
        else:
            return self.coordinates[self.current]

    def __getitem__(self, key):
        if key >= self.dimension:
            raise IndexError
        else:
            return self.coordinates[key]

    def __str__(self):
        return 'Vector: {}'.format(self.coordinates)

    def __eq__(self, v, tolerance=1e-10):
        result = True
        for i in range(self.dimension):
            result = (
                result & (abs(self.coordinates[i] - v.coordinates[i]) < tolerance))
        return result

    def plus(self, v):
        new_coordinates = [x + y for x,
                           y in zip(self.coordinates, v.coordinates)]
        return Vector(new_coordinates)

    def minus(self, v):
        new_coordinates = [x - y for (x, y) in zip(self.coordinates, v.coordinates)]
        return Vector(new_coordinates)

    def scalar_multiply(self, c):
        new_coordinates = [Decimal(c) * x for x in self.coordinates]
        return Vector(new_coordinates)

    def magnitude(self):
        coordinates_squared = [x**Decimal(2) for x in self.coordinates]
        return Decimal(math.sqrt(sum(coordinates_squared)))
        # math.sqrt的结果是float啊

    def normalization(self):  # make a unit vector
        try:
            magnitude = self.magnitude()
            return self.scalar_multiply(Decimal(1.0) / magnitude)
        except ZeroDivisionError:
            raise Exception(self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG)

    def dot_product(self, v):
        return sum([x * y for (x, y) in zip(self.coordinates, v.coordinates)])

    def angel_of_two_vectors(self, v, in_degrees=False):
        try:
            u1 = self.normalization()  # unit vector
            u2 = v.normalization()
            angel_in_radians = math.acos(round(u1.dot_product(u2),5))
            # 注意这里dot的结果可能大于1或小于-1，需要四舍五入
            if in_degrees:
                degrees_per_radian = 180. / math.pi
                return angel_in_radians * degrees_per_radian
            else:
                return angel_in_radians

        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception('Cannot compute an angle with the zero vector')
            else:
                raise e

    def parallelism(self, v, tolerance=1e-10):
        return (self.magnitude() < tolerance or
                v.magnitude() < tolerance or
                self.angel_of_two_vectors(v) == 0 or
                self.angel_of_two_vectors(v) == math.pi)

    def orthogonality(self, v, tolerance=1e-10):
        return abs(self.dot_product(v)) < tolerance

    def projection_vector(self, b):
        try:
            ub = b.normalization()
            weight = self.dot_product(ub)
            return ub.scalar_multiply(weight)
        except Exception as e:  # when basis vector is 0
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception(self.NO_UNIQUE_PARALLEL_COMPONENT_MEG)
            else:
                raise e

    def orthogonal_vector(self, b):
        try:
            projection = self.projection_vector(b)
            return self.minus(projection)
        except Exception as e:
            if str(e) == self.NO_UNIQUE_PARALLEL_COMPONENT_MEG:
                raise Exception(self.NO_UNIQUE_PARALLEL_COMPONENT_MEG)
            else:
                e

    def cross_product(self, v):
        try:
            x_1, y_1, z_1 = self.coordinates
            x_2, y_2, z_2 = v.coordinates
            new_coordinates = [y_1 * z_2 - y_2 * z_1,
                               -(x_1 * z_2 - x_2 * z_1),
                               x_1 * y_2 - x_2 * y_1]
            return Vector(new_coordinates)
        except ValueError as e:
            msg = str(e)
            if msg == 'need more than 2 values to unpack':
                self_embedded_in_R3 = Vector(self.coordinates + ('0',))
                v_embedded_in_R3 = Vector(v.coordinates + ('0',))
                return self_embedded_in_R3.cross_product(v_embedded_in_R3)
            elif (msg == 'too many values to unpack' or
                  msg == 'need more than 1 value to unpack'):
                raise Exception(self.ONLY_DEFINE_IN_TWO_THREE_DIMS_MSG)
            else:
                raise e

    def area_of_parallelogram(self, v):
        return self.cross_product(v).magnitude()

    def area_of_triangle(self, v):
        return self.area_of_parallelogram(v) / Decimal(2.)

    def is_zero(self, tolerance=1e-10):
        result = True
        for x in self.coordinates:
            if abs(x) < tolerance:
                result &= True
            else:
                result &= False
                break
        return result
