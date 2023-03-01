import sympy as sp


def calc_miller_indices(intercepts):
    """Calculate Miller indices for a set of intercepts.

    Parameters
    ----------
    intercepts : list
        The intercepts of the crystal planes.

    Returns
    -------
    miller_indices : list
        The Miller indices of the crystal planes.

    """
    #replace any zeros with 0.000000000001
    intercepts = [intercept if intercept != 0 else 0.000000000001 for intercept in intercepts]
    reciprocals = [1/intercept for intercept in intercepts]
    #make smallest as 1
    smallest = min(reciprocals)
    miller_indices = [reciprocal/smallest for reciprocal in reciprocals]
    #check if they are integers
    if all([miller_index.is_integer() for miller_index in miller_indices]):
        idxs = [int(x) for x in miller_indices]
        #replace any number grater than 1000 with 'inf'
        idxs = [idx if idx < 1000 else 'inf' for idx in idxs]
        print(f'Miller indices: {idxs}')
        return idxs
    else:
        #find what the smallest non-integer is
        smallest_non_integer = min([miller_index for miller_index in miller_indices if 
            not miller_index.is_integer()])
        smallest_non_integer_fraction = sp.Rational(smallest_non_integer).limit_denominator()
        #multiply all miller indices by the denominator
        miller_indices = [miller_index*smallest_non_integer_fraction.denominator for
            miller_index in miller_indices]
        idxs = [int(x) for x in miller_indices]
        #replace any number grater than 1000 with 'inf'
        idxs = [idx if idx < 1000 else 'inf' for idx in idxs]
        print(f'Miller indices: {idxs}')
        return idxs


