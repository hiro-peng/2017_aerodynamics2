# -*- coding: utf-8 -*-
"""
Created on Fri Jul 28 17:01:34 2017
@author: folding
"""
import matplotlib.pyplot as plt
import math
import cmath
import numpy as n

"""翼型パラメータ"""
a = 0.6 #写像化する前の円の半径→翼の厚さに相当
v0 = 1.0 #風速
deg = math.pi/180.0
alfa = 5.0*deg #迎角
beta = 20.0*deg #キャンバーの反りの角度
c = 0.5 #ずらす距離
ganma = 4.0*math.pi*v0*a*math.sin(alfa + beta) #クッタの条件

"""翼のパラメータを複素数に"""
cv = complex(v0,0.0) #流速
ca = complex(a,0.0) #写像前の円の半径
cganma = complex(0.0,ganma/(2.0*math.pi)) #ix循環/2π
Zc = complex(c,0.0) + ca*cmath.exp(complex(0.0,math.pi-beta)) #Zcを求める

"""プロットパラメータ"""
max_theta = 180
max_r = 2500
dtheta = 2.0*deg
dr = 0.001

"""格納する配列作成"""
x = np.ones((max_theta,max_r))
y = np.ones((max_theta,max_r))
f = np.zeros((max_theta,max_r))
g = np.zeros((max_theta,max_r))
xi = np.zeros((max_theta,max_r))
eta = np.zeros((max_theta,max_r))

def polar():
    for theta in range(max_theta):
        for r in range(max_r):

            """座標設定"""
            r1 = a + dr*r
            theta1 = -beta + dtheta*theta

            x[theta][r] = r1*math.cos(theta1)
            y[theta][r] = r1*math.sin(theta1)

            Z = complex(x[theta][r],y[theta][r])
            Z0 = Z + Zc

            """以下、座標を代入"""
            #小文字zを求める↓
            cz = Z *cmath.exp(complex(0.0,-alfa))

            #半径aで循環Γを持つポテンシャル流れ関数f↓
            cf = cv*(cz + ca*ca/cz) + cganma*cmath.log(cz)

            #流線関数
            f[theta][r] = cf.real
            g[theta][r] = cf.imag

            """翼型決定"""
            zeta = Z0 + c*c/Z0
            xi[theta][r] = zeta.real
            eta[theta][r] = zeta.imag

def stream_plot():
    """範囲指定"""
    plt.xlim(-2.0,2.0)
    plt.ylim(-2.0,2.0)

    """流線を描く"""
    plt.contour(xi[:,:],eta[:,:],g[:,:],locator = plt.MultipleLocator(0.05))
    plt.summer()

    """翼を描く(翼はξ、ηの配列のr=0の部分)"""
    plt.plot(xi[:,0],eta[:,0])

    plt.title('Sream Line')
    plt.xlabel("ξ")
    plt.ylabel("η")
    cbar = plt.colorbar()
    cbar.set_label("ψ")
    plt.show()


if __name__ == "__main__":
    polar()
    stream_plot()
