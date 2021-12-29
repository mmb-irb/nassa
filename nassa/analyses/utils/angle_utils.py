def fix_angle_range(x, domain=[0, 360]):
    """Fix angle range so it's in the given angle range (degrees) by adding or subtracting 360.

    :param float x: angle value (asumed to be in degrees)
    :param sequence domain: start and end of angle range.

    : return float: angle value with fixed range """
    while x < domain[0]:
        x += 360
    while x > domain[1]:
        x -= 360
    return x
