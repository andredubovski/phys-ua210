from math import sqrt

class Cell:

    #constants
    gamma = 1.4

    #initial conditions
    rho = 0
    v = 0
    e = 0
    E = 0
    P = 0
    lambda_plus = 0
    lambda_minus = 0

    #variables
    alpha_plus = [0, 0, 0]
    alpha_minus = [0, 0, 0]

    def __init__(self, rho, v, e) -> None:
        self.rho = rho
        self.v = v
        self.e = e
        self.P = (self.gamma-1)*self.rho*self.e
        self.E = rho*e + 1/2*rho*v**2.0
        self.U = [rho, rho*v, self.E]
        self.F = [rho*v,
                  rho*(v**2.0+(self.gamma - 1)*e),
                  rho*v*(1/2*v**2+self.gamma*e)]
        print(self.F)
        self.lambda_plus = v + sqrt((self.gamma*(self.gamma-1)*self.rho*self.e)/self.rho)
        self.lambda_minus = v - sqrt((self.gamma*(self.gamma-1)*self.rho*self.e)/self.rho)
        print("v:", v)
        print("lambda_plus:", self.lambda_plus)
        print("lambda_minus:", self.lambda_minus)

    
    def find_F_half(self, left_cell):

        for i in range(3):
            left_cell.alpha_plus[i] = max([0, left_cell.lambda_plus*self.U[i], left_cell.lambda_plus*self.U[i]])
            left_cell.alpha_minus[i] = max([0, -left_cell.lambda_minus*self.U[i], -left_cell.lambda_minus*self.U[i]])

        for i in range(3):
            left_cell.F_half[i] = (left_cell.alpha_plus[i]*left_cell.F[i] + left_cell.alpha_minus[i]*self.F[i] \
            - left_cell.alpha_plus[i]*left_cell.alpha_minus[i]*(self.U[i]-left_cell.U[i])) \
            / (left_cell.alpha_plus[i] + left_cell.alpha_minus[i])

        return self.F_half
    
    def find_F(self):
        self.P = (self.gamma-1)*self.U[0]*(self.U[2]/self.U[0]-1/2*(self.U[1]/self.U[0])**2)
        self.F[0] = self.U[1]
        self.F[1] = self.U[1]**2/self.U[0]+self.P
        self.F[2] = (self.U[2]+self.P)*self.U[1]/self.U[0]
        
        return self.F
    
    def evolve(self, delta_t, delta_x, left_cell):
        for i in range(3):
            self.U[i] = self.U[i] - delta_t*(self.F_half[i] - left_cell.F_half[i])/delta_x
        print(self.F_half[0])
        print(left_cell.F_half[0])
        print('-')
        self.rho = self.U[0]
        self.v = self.U[1]
        self.F = self.find_F()
        


    U = [0, 0, 0]
    F = [0, 0, 0]
    F_half = [0, 0, 0]

    def update():
        return