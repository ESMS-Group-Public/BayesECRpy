import numpy as np

def Rejection(variable_min, variable_max, variable, SIGMA_variable):

    LogOf_variable_min = np.log10(variable_min)
    LogOf_variable_max = np.log10(variable_max)
    variable_Flag = True

    LogOf_variable_0 = np.log10(variable)

    while variable_Flag == True:
        LogOf_variable_1pre = LogOf_variable_min + (LogOf_variable_max - LogOf_variable_min)* np.random.uniform(0,1,(1,1))
        N_variable = -1*(LogOf_variable_1pre - LogOf_variable_0)**2
        P_variable = np.exp(N_variable / (2 * SIGMA_variable ** 2)) / (SIGMA_variable * (2 * np.pi) ** .5)
        Pr0_variable = 1 / (SIGMA_variable * (2 * np.pi) ** .5)
        PROB_variable = P_variable / Pr0_variable
        if (PROB_variable > np.random.uniform()):
            LogOf_variable_0 = LogOf_variable_1pre
            variable_Flag = False

    proposed_variable = 10**LogOf_variable_0

    return proposed_variable


