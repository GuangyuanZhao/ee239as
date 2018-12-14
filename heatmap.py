import matplotlib.pyplot as plt
import numpy as np
import glob as glob
import cv2
import random

# Hyperparameter
# Window size
WINDOW_WIDTH = 400
WINDOW_HEIGHT = 400


# Define polar coordinates
def cart2pol(y, x, center=(0,0)):
    rho = np.sqrt((y - center[0]) ** 2 + (x - center[1]) ** 2)
    phi = np.arctan2(y - center[0], x - center[1])
    if phi < 0:
        phi += 2*np.pi
    return rho, phi

def pol2cart(rho, phi, center=(0,0)):
    y = round(rho * np.sin(phi) + center[0])
    x = round(rho * np.cos(phi) + center[1])
    return int(y), int(x)

class Aop_Map:
    def __init__(self, prefix):
        """
        :param prefix: prefix should not end up with /
        """
        self.prefix = prefix
        self.gray_images = None
        self.aop = None
        self.dolp = None
        self.zc = None

    def read_4angles(self):
        self.gray_images = []
        for angle in [0,45,90,135]:
            dir = self.prefix + '/%d/*.JPG'%angle
            files = glob.glob(dir)
            img = []
            for i in range(len(files)):
                print("Reading image %d for angle %d" % (i+1, angle))
                single_img = plt.imread(files[i])
                single_img = np.asarray(0.299 * single_img[..., 0] + 0.587 * single_img[..., 1] + 0.114 * single_img[..., 2])
                img.append(single_img)
            self.gray_images.append(np.average(img, axis=0))

        return self.gray_images

    def compute_dolp(self):
        if self.gray_images is None:
            self.read_4angles()
        s0 = np.sum(self.gray_images, axis=0)/4.0
        s1 = self.gray_images[0] - self.gray_images[2]
        s2 = self.gray_images[1] - self.gray_images[3]

        s12 = [s1, s2]
        self.dolp = np.linalg.norm(s12, axis=0)/s0
        return self.dolp

    def compute_aop(self):
        if self.gray_images is None:
            self.read_4angles()
        # Compute stokes vector
        s0 = np.sum(self.gray_images, axis=0)/4.0
        s1 = self.gray_images[0] - self.gray_images[2]
        s2 = self.gray_images[1] - self.gray_images[3]

        # Compute pixelwise aop
        self.aop = 0.5*np.arctan2(abs(s2), abs(s1))
        return self.aop

    def get_zenith_center(self):
        if self.aop is None:
            self.compute_aop()
        self.grad = cv2.Laplacian(self.aop, cv2.CV_64F)
        coor = np.where(self.grad==self.grad.max())
        coor_y = round(np.average(coor[0]))
        coor_x = round(np.average(coor[1]))
        self.zc = (int(coor_y), int(coor_x))
        return self.zc

    def get_zenith_center_alt(self, y1=0, x1=0, y2=3000, x2=4000):
        if self.dolp is None:
            self.compute_dolp()
        coor = np.where(self.dolp[y1:y2, x1:x2]==self.dolp[y1:y2, x1:x2].min())
        coor_y = round(np.average(coor[0]) + y1)
        coor_x = round(np.average(coor[1]) + x1)
        self.zc = (int(coor_y), int(coor_x))
        return self.zc


    def get_random_coor(self, aop, num, th=0.5):
        h, w = aop.shape
        x_coor = []
        y_coor = []
        for n in range(num):
            y = random.randint(0, h-1)
            x = random.randint(0, w-1)
            while aop[y, x] < th:
                x = random.randint(0, w-1)
            y_coor.append(y)
            x_coor.append(x)
        return np.asarray(y_coor), np.asarray(x_coor)


    def compute_azimuth_line(self, center=None):
        # get zenith center and Aop
        if center is None:
            if self.zc is None:
                center = self.get_zenith_center()
            else:
                center = self.zc
        if self.aop is None:
            self.compute_aop()

        # Get pixel indices on the half circle
        def _get_circle(self, r=100):
            half_circle = {}
            zc = np.asarray(self.zc, dtype=int)
            x_range = list(range(zc[1]-r, zc[1]+r+1))
            for x in x_range:
                y = zc[0] + np.sqrt(r**2 - (x-zc[1])**2)
                half_circle[x] = int(y)
            return half_circle

        radius = int(min(WINDOW_WIDTH, WINDOW_HEIGHT)/2)
        clc_dict = _get_circle(self, radius)


        # Random sampling 1000 points
        x_dim_max = int(WINDOW_WIDTH/2)
        y_dim_max = int(WINDOW_HEIGHT/2)

        y_list, x_list = self.get_random_coor(aop[center[0]-y_dim_max: center[0]+y_dim_max+1,
                             center[1]-x_dim_max: center[1]+x_dim_max+1], num=1000, th=0.5)
        y_list = y_list + center[0]-y_dim_max
        x_list = x_list + center[1]-x_dim_max


        # Optimize angle by angle
        opt_sol = (0,0)
        opt_err = None
        for x0 in sorted(clc_dict.keys(), reverse=True):
            y0 = clc_dict[x0]
            if self.aop[y0, x0] < 0.5:
                continue
            _, phi0 = cart2pol(y0, x0, center)
            error = 0
            for i in range(len(x_list)):
                # Compute the symmetric index
                rou, phi = cart2pol(y_list[i], x_list[i], center)
                y_sym, x_sym = pol2cart(rou, 2*phi0-phi, center)
                error += abs(aop[y_list[i], x_list[i]]-aop[y_sym, x_sym])

            if opt_err is None:
                opt_err = error
                opt_sol = (y0, x0)
            else:
                if error < opt_err:
                    opt_err = error
                    opt_sol = (y0, x0)
                elif error == opt_err:
                    print("Multiple opt found! Error: %f, pixel: (%d, %d)" % (error, y0, x0))

        self.az_line = opt_sol

        # Draw the line
        p1 = [int(2*opt_sol[1]-1*center[1]), int(2*opt_sol[0]-1*center[0])]
        p2 = [int(3*center[1]-2*opt_sol[1]), int(3*center[0]-2*opt_sol[0])]

        self.az_pic = np.copy(self.aop)
        self.az_pic = cv2.line(self.az_pic, (p1[0], p1[1]), (p2[0], p2[1]), color = (0,255,0), thickness=20)
        plt.imshow(self.az_pic)

        return self.az_line


if __name__ == '__main__':
    # map1 = Aop_Map(prefix='./11_31/group1')
    # center = map1.get_zenith_center_alt(1414, 1762, 2486, 2850)

    
    # map1 = Aop_Map(prefix='./12_01/group5')
    map1 = Aop_Map(prefix='./12_13/13_20')
    aop = map1.compute_aop()

    center = map1.get_zenith_center_alt()
    az_line = map1.compute_azimuth_line(center)

    theta1 = np.rad2deg(cart2pol(az_line[0], az_line[1], center)[1])






