# coding=utf-8

import math
from decimal import Decimal, getcontext

getcontext().prec = 15


class Vector(object):

    CANNOT_NORMALIZE_ZERO_VECTOR_MSG = 'Cannot normalize the zero vector'

    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple([Decimal(x) for x in coordinates])
            self.dimension = len(coordinates)

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')

    def __str__(self):
        return 'Vector: {}'.format(self.coordinates)

    def __eq__(self, v):
        return self.coordinates == v.coordinates

    def plus(self, v):
        new_coordinates = [x + y for x,
                           y in zip(self.coordinates, v.coordinates)]
        # print zip(self.coordinates, v.coordinates)
        #[(8.218, -1.129), (-9.341, 2.111)]
        return Vector(new_coordinates)

    def minus(self, v):
        new_coordinates = [
            x - y for (x, y) in zip(self.coordinates, v.coordinates)]
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
            angel_in_radians = math.acos(u1.dot_product(u2))

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

# my own code
    def parallelism(self, v, tolerance=1e-10):
        return (self.magnitude() < tolerance or
                v.magnitude() < tolerance or
                self.angel_of_two_vectors(v) == 0 or
                self.angel_of_two_vectors(v) == math.pi)

    def orthogonality(self, v, tolerance=1e-10):
        return abs(self.dot_product(v)) < tolerance

# test for parallel and orthogonal vectors
# v and w are parallel if one is a scalar multiply of the other
# v and w are orthogonal if their dot product is zero
# 0:parallel and orthogonal to all vectors
my_vector11 = Vector([-7.579, -7.88])
my_vector12 = Vector([22.737, 23.64])
my_vector21 = Vector([-2.029, 9.97, 4.172])
my_vector22 = Vector([-9.231, -6.639, -7.245])
my_vector31 = Vector([-2.328, -7.284, -1.214])
my_vector32 = Vector([-1.821, 1.072, -2.94])
my_vector41 = Vector([2.118, 4.827])
my_vector42 = Vector([0, 0])

print my_vector11.parallelism(my_vector12)
print my_vector11.orthogonality(my_vector12)


print my_vector21.parallelism(my_vector22)
print my_vector21.orthogonality(my_vector22)

print my_vector31.parallelism(my_vector32)
print my_vector31.orthogonality(my_vector32)

print my_vector41.parallelism(my_vector42)
print my_vector41.orthogonality(my_vector42)
