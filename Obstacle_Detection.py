from math import atan, sqrt, pi, pow
from pyfirmata import Arduino, SERVO
from time import sleep


def detect_obstacle(dir_vec):
    """ This method uses an ultrasonic range sensor
        and servo motor to detect obstacles in the path
        of an individual based on their direction of travel.
    Args:
        dir_vec (list): A list containing the x, y, and z coordinates
                        retrieved from the accelerometer data.
    Returns:
        Nothing. A piezoelectric buzzer is triggered if an obstacle is detected
        within a threshold of 50 cm.
        Note: For now, the yaw angle is being returned.
    """

    ax = dir_vec[0]
    ay = dir_vec[1]
    az = dir_vec[2]

    pitch = 180 * atan(ax / sqrt(pow(ay, 2) + pow(az, 2))) / pi
    roll = 180 * atan(ay / sqrt(pow(ax, 2) + pow(az, 2))) / pi

    # We want to calculate the yaw angle to determine the servo rotation.
    # Using pythagoras theorem and converting from radians to degrees.
    yaw = 180 * atan(az / sqrt(pow(ax, 2) + pow(az, 2))) / pi

    return(yaw)

    # TODO: Experiment with code below to see if controlling the ultrasonic sensor
    # from python is possible (so far there are some brick walls)
    '''
    # Setup Arduino UNO for obstacle detection.
    port = 'COM3'
    board = Arduino(port)
    # Wait to let the Arduino board and pyfirmata synchronize.
    sleep(5)

    # Run once
    loop = 1
    for x in range(int(loop)):
        # Assign pins to variables.
        trig = board.get_pin('d:10:o')
        echo = board.get_pin('d:11:i')

        # Use before reading.
        echo.enable_reporting()

        # Rotate the servo to the angle determined at the start of this method.
        board.digital[12].mode = SERVO
        board.digital[12].write(yaw)
        sleep(0.015)


        # TODO: Reset the servo to original position after rotating and measuring distance!
    '''

