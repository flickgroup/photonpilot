import numpy as np

class LayeredDGF:
    """
    Dyadic Green's Tensor class for a layered medium
    This implementation only works at r=r'
    """
    def __init__(self, qs, k, r12_s, r13_s, r12_p, r13_p, d, b, t):
        """
        Init function for the LayeredDGF class
        :param qs: In-plane wave number array for integration
        :param k: wave number inside an infinite medium made of the 
        same material as the cavity layer
        :param r12_s: Fresnel coefficient of the bottom material-mirror 
        interface for s-polarized light
        :param r13_s: Fresnel coefficient of the top material-mirror 
        interface for s-polarized light
        :param r12_p: Fresnel coefficient of the bottom material-mirror 
        interface for p-polarized light
        :param r13_p: Fresnel coefficient of the top material-mirror 
        interface for p-polarized light
        :param d: thickness of the cavity region
        :param b: distance from the emitter to the bottom mirror
        :param t: distance from the emitter to the top mirror
        """
        self.qs = qs
        self.k = k
        self.kz = np.sqrt(self.k**2 - self.qs**2 + 1E-3)
        self.dq = np.abs(self.qs[1] - self.qs[0])
        self.r12_s = r12_s
        self.r13_s = r13_s
        self.r12_p = r12_p
        self.r13_p = r13_p
        self.d = d
        self.b = b
        self.t = t
        self.get_Fzz_te_ip()
        self.get_Fzz_tm_ip()
        self.get_Fzz_tm_op()

    def get_all_DGF_components(self):
        """
        Helper function to get all DGF components
        """
        self.get_DGF_zz_te_ip()
        self.get_DGF_zz_tm_ip()
        self.get_DGF_zz_tm_op()

    def get_DGF_zz_te_ip(self):
        """
        Function to calculate the z = z' component of the
        DGF in the layered cavity setup
        """
        term1 = 1j/(8*np.pi)
        term2 = self.qs/self.kz*self.dq
        self.G_zz_te_ip = term1*np.sum(term2*self.Fzz_te_ip)

    def get_DGF_zz_tm_ip(self):
        """
        Function to calculate the z = z' component of the
        DGF in the layered cavity setup
        """
        term1 = 1j/(self.k**2*8*np.pi)
        term2 = self.qs*self.kz*self.dq
        self.G_zz_tm_ip = term1*np.sum(term2*self.Fzz_tm_ip)

    def get_DGF_zz_tm_op(self):
        """
        Function to calculate the z = z' component of the
        DGF in the layered cavity setup
        """
        term1 = 1j/(self.k**2*8*np.pi)
        term2 = 2*self.qs**3/self.kz*self.dq
        self.G_zz_tm_op = term1*np.sum(term2*self.Fzz_tm_op)

    def get_Fzz_te_ip(self):
        """
        Function to get the out-of-plane modulation
        functions at z = z' for horizontal dipole
        orientation and TE-polarized light
        :return Fzz_sh: out-of-plane modulation function
        for the layered cavity
        """
        term1 = (1 + self.r12_s*np.exp(2j*self.kz*self.b))
        term2 = (1 + self.r13_s*np.exp(2j*self.kz*self.t))
        term3 = (1 - self.r12_s*self.r13_s*np.exp(2j*self.kz*self.d) + 0.01)
        self.Fzz_te_ip = term1*term2/term3

    def get_Fzz_tm_ip(self):
        """
        Function to get the out-of-plane modulation
        functions at z = z' for horizontal dipole
        orientation and TM-polarized light
        :return Fzz_ph: out-of-plane modulation function
        for the layered cavity
        """
        term1 = (1 - self.r12_p*np.exp(2j*self.kz*self.b))
        term2 = (1 - self.r13_p*np.exp(2j*self.kz*self.t))
        term3 = (1 - self.r12_p*self.r13_p*np.exp(2j*self.kz*self.d) + 0.01)
        self.Fzz_tm_ip =  term1*term2/term3

    def get_Fzz_tm_op(self):
        """
        Function to get the out-of-plane modulation
        functions at z = z' for vertical dipole
        orientation and TM-polarized light
        :return Fzz_ph: out-of-plane modulation function
        for the layered cavity
        """
        term1 = (1 + self.r12_p*np.exp(2j*self.kz*self.b))
        term2 = (1 + self.r13_p*np.exp(2j*self.kz*self.t))
        term3 = (1 - self.r12_p*self.r13_p*np.exp(2j*self.kz*self.d) + 0.01)
        self.Fzz_tm_op =  term1*term2/term3


