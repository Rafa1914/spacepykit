class TLE:
    def __init__(
        self, a: float, i: float, ra: float, e: float, w: float, theta: float
    ) -> None:
        self.a = a
        self.i = i
        self.ra = ra
        self.e = e
        self.w = w
        self.theta = theta
    
    def to_list(self):
        return [self.a,self.i,self.ra,self.e,self.w]