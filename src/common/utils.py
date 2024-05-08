def convert_angle_0_to_360_range(theta:float) -> float:
    '''
    Convert the angle theta in a [0,2pi] angle range.

    Attributes:
    theta: angle to be converted [deg]

    Returns:
    Converted angle [deg]
    '''
    if 0 <= theta <= 360:
        theta = theta
    elif theta > 360:
        n = int(theta/360) #Gets the integer part of the quocient
        theta = theta - 360*n
    else:
        n = abs(int(theta/360))
        theta =  (n+1)*360 + theta
    return theta