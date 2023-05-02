import numpy as np

class SphericalDGF:
    """
    Class to hold the Dyadic Green's tensor of the spherical 
    cavity at r=r'

    The current implementation only works for emitter positions
    at the center of the cavity such that r=r'=0 
    """
    def __init__(self, R, omega, epsilon_shell):
        """
        Init function for the SphericalDGF class
        :param R: Inner radius of the spherical cavity
        :param omega: Array of frequencies for which to
        calculate the cavity field strengths
        :param epsilon3: frequency dependent dielectric
        function of the cavity region'
        :param epsilon2: frequency dependent dielectric
        function of the mirror region
        """
        # Fundamental constants
        self.c = 3E8
        self.hbar = 1.05E-34
        self.e = 1.6E-19
        self.epsilon0 = 8.9E-12

        # Setting inputs
        self.R = R
        self.omega = omega
        self.dw = np.abs(np.diff(omega)[0])
        self.epsilon_shell = epsilon_shell
        self.n_shell = np.sqrt(self.epsilon_shell)
        
        # Calculating size paramters
        self.rho=self.R*self.omega/self.c

    def get_lambda2(self):
        """
        Function to get the squared norm of the cavity field 
        strengths at r=r'=0 
        """
        e = self.e
        dw = self.dw
        omega = self.omega
        epsilon0 = self.epsilon0
        c = self.c

        self.get_reflection_coefficient()
        s = 1 + np.real(self.r22p)
        self.lambda2 = 2*dw*e**2*omega**2/(6*np.pi**2*epsilon0*c**3)*s
        
        # Converting from J/m² to to eV/m²
        self.lambda2 *= (6.24E18)/1E18

    def get_DGF(self):
        """
        Helper function to get the DGF
        """
        self.get_reflection_coefficient()
        self.DGF = 1j*self.omega/(6*np.pi*self.c)*(1 + self.r22p)

    def get_reflection_coefficient(self):
        """
        Function to get the reflection coefficient for the 
        lowest order TM mode of the cavity
        """
        rho = self.rho
        n = self.n_shell

        # Calculating the reflection coefficient according to SI Eq. B3
        numerator = np.exp(1j*rho)*(1j + rho*(n +1) - 1j*rho**2*n - rho**3*n**2/(n+1))
        denominator = np.sin(rho) - rho*(np.cos(rho) + 1j*n*np.sin(rho)) + 1j*rho**2*n*np.cos(rho) - rho**3*(np.cos(rho) - 1j*n*np.sin(rho))*n**2/(n**2-1)

        self.r22p = numerator/denominator


