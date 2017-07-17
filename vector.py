# coding=utf-8

import math
from decimal import Decimal, getcontext

getcontext().prec = 8


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
            return self.scalar_multiply(Decimal(1) / magnitude)
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

# def plus(v1, v2):
#     result = []
#     for i in range(v1.dimension):
#         result.append(v1.coordinates[i] + v2.coordinates[i])
#         i += 1
#     return result


# def minus(v1, v2):
#     result = []
#     for i in range(v1.dimension):
#         result.append(v1.coordinates[i] - v2.coordinates[i])
#         i += 1
#     return result


# def scalar_multiply(v1, n):
#     result = []
#     for i in range(v1.dimension):
#         result.append(v1.coordinates[i] * n)
#         i += 1
#     return result
    # def magnitude(self):
    #     result = 0
    #     for i in self.coordinates:
    #         result += i * i
    #     result = math.sqrt(result)
    #     return result

    # def direction(self):
    #     new_coordinates = [(1/self.magnitude()) * x for x in self.coordinates]
    #     return Vector(new_coordinates)
    # def dot_product(self, v):
    #     return sum(x * y for (x, y) in zip(self.coordinates, v.coordinates))

    # def angel_of_two_vectors(self, v):
    # return math.acos(self.dot_product(v) / self.magnitude() / v.magnitude())
    def parallelism(self, v):
        
        if sum([x==0 for x in v.coordinates])==True:
            return True
        else:
            divide = self.coordinates[0]/v.coordinates[0]
            if self.scalar_multiply(divide[0]).coordinates == self.coordinates:
                return True
            else:
                return False

    def orthogonality(self, v):
        if self.dot_product(v) == 0:
            return True
        else:
            return False


# test for plus,minus,scalar_multiply,__str__(),__eq__()
# my_vector1 = Vector([8.218, -9.341])
# my_vector2 = Vector([-1.129, 2.111])
# my_vector3 = Vector([7.119, 8.215])
# my_vector4 = Vector([-8.223, 0.878])
# my_vector5 = Vector([1.671, -1.012, -0.318])
# scalar = 7.41
# # print my_vector1.__str__()
# # print my_vector1.__eq__(my_vector2)

# print my_vector1.plus(my_vector2)
# print my_vector3.minus(my_vector4)
# print my_vector5.scalar_multiply(scalar)


# # test for magnitude and direction

# test for inner products(dot) & angel
# my_vector11 = Vector([7.887, 4.138])
# my_vector12 = Vector([-8.802, 6.776])
# my_vector21 = Vector([-5.955, -4.904, -1.874])
# my_vector22 = Vector([-4.496, -8.755, 7.103])
# my_vector31 = Vector([3.183, -7.627])
# my_vector32 = Vector([-2.668, 5.319])
# my_vector41 = Vector([7.35, 0.221, 5.188])
# my_vector42 = Vector([2.751, 8.259, 3.985])

# print my_vector11.dot_product(my_vector12)
# print my_vector21.dot_product(my_vector22)
# print my_vector11.magnitude()
# print my_vector11.normalization()
# print my_vector31.angel_of_two_vectors(my_vector32)
# print my_vector41.angel_of_two_vectors(my_vector42, in_degrees=True)
# print math.pi

# test for parallel and orthogonal vectors
# v and w are parallel if one is a scalar multiply of the other
# v and w are orthogonal if their dot product is zero
# 0:parallel and orthogonal to all vectors
my_vector11 = Vector([-7.579, -7.88])
my_vector12 = Vector([22.737, 23.64])
my_vector21 = Vector([-2.029, 9.97, 4.172])
my_vector22 = Vector([-9.231, -6.639, -7.245])
my_vector31 = Vector([-2.328, -7.284, -1.24])
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