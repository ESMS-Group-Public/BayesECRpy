import numpy as np

def METROPOLIS(proposedLikely, LOGLIKELY, fact):
    '''
    METROPOLIS accepts or rejects a drawn candidate

    :return: outputFlag
    '''

    someNum = fact*(proposedLikely - LOGLIKELY)
    P = np.log(np.random.rand())
    if someNum > 0:
        outputFlag = 1
    elif someNum > P:
        outputFlag = 1
    else:
        outputFlag = 0

    return outputFlag

