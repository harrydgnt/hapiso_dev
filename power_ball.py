import numpy as np
import sys
from itertools import combinations, permutations
def number():
    number_array=list(range(1,76))
    rand_num = np.random.randint(1,100000000000000)
    marker=0
    index=rand_num%11238513
    white_balls=[]
    for item in combinations(number_array, 5):
        marker = marker +1
        if marker == index:
            white_balls=item
    power_ball = np.random.randint(1,26)
    return white_balls, power_ball

def run():
    winning_set, winning_power_ball = number()


if __name__ == "__main__":
    print "hello world "
    run()