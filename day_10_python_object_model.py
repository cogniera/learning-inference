import math

class HyperboloidPoint:
    def __init__(self, t, x=0, y=0):
        # inputs for the point
        self.t = t 
        self.x = x
        self.y = y
        self.origin = (1, 0, 0)    

        # Calculating the minkowski self product
        minkowski_self_product = self._calculate_minkowski_self_product()
                
        # Validation logic
        if not math.isclose(minkowski_self_product, -1.0, rel_tol=1e-9):
            raise ValueError("Point is not on the hyperboloid")
            
        # the time should be forward
        if self.t <= 0:
            raise ValueError("Time coordinate must be positive") 
    
    def _calculate_minkowski_self_product(self):
        return (math.pow(self.x, 2) + math.pow(self.y, 2) - math.pow(self.t, 2))
    
    def _calculate_minkowski_inner_product(self, other):
        return (-self.t * other.t + self.x * other.x + self.y * other.y)
    
    def __repr__(self):
        # print the coordinates in the same order as __init__
        return f"HyperboloidPoint({self.t:.3f}, {self.x:.3f}, {self.y:.3f})"
    
    def __matmul__(self, other):
        # matmul is the minkowski inner product
        minkowski_inner_product = self._calculate_minkowski_inner_product(other)
        
        return minkowski_inner_product
    
    def hyperbolicDistance(self, other):
        # calculate the inner product
        minkowski_inner_product = self._calculate_minkowski_inner_product(other)
        
        # distance is just arccosh of this inner product
        dist = math.acosh(-minkowski_inner_product)
        
        return dist
    
    def __abs__(self):
        return math.acosh(self.t)
    