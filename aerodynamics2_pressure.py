# -*- coding: utf-8 -*-
"""
Created on Sun Jul 30 16:29:19 2017

@author: folding
"""
import matplotlib.pyplot as plt
import math
import cmath
import numpy as np

"""翼型パラメータ"""
a = 0.6 #写像化する前の円の半径→翼の厚さに相当
v0 = 1.0 #風速
roh = 1.0 #空気密度
deg = math.pi/180.0
beta = 20.0*deg #キャンバーの反りの角度
c = 0.5
l = 4.0*c #ずらす距離(コード長の１/４)
alfa = 5.0*deg
5.0*deg

"""翼のパラメータを複素数に"""
cc = complex(c,0.0) #ずらす距離(コード長の１/４)
cv = complex(v0,0.0) #風速
ca = complex(a,0.0) #写像前の円の半径
Zc = complex(c,0.0) + ca*cmath.exp(complex(0.0,math.pi-beta)) #Zcを求める

"""プロットパラメータ"""
max_theta = 180
max_r = 10**4
dtheta = 2.0*deg
dr = 0.001

"""格納する配列作成"""
x = np.ones((max_theta,max_r))
y = np.ones((max_theta,max_r))
f = np.zeros((max_theta,max_r))
g = np.zeros((max_theta,max_r))
xi = np.zeros((max_theta,max_r))
eta = np.zeros((max_theta,max_r))
Cp0 = np.zeros((max_theta,max_r))

"""メッシュ、変換"""
def polar():

    """迎角αによって変化する条件"""
    ganma = 4.0*math.pi*v0*a*math.sin(alfa + beta) #クッタの条件
    cganma = complex(0.0,ganma/(2.0*math.pi)) #ix循環/2π

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
            #ポテンシャル流れ関数fのζでの微分↓
            dcf = cmath.exp(complex(0.0,-alfa))*(cv*(complex(1.0,0.0) - (ca*ca) / (cz*cz)) + (cganma/cz)) / (complex(1.0,0.0) - (cc*cc)/(Z0*Z0))
            #圧力無次元化↓
            Cp0[theta][r] = 1 - (((dcf.real**2)+(dcf.imag**2)) / (v0**2))

            """翼型決定"""
            zeta = Z0 + c*c/Z0
            xi[theta][r] = zeta.real
            eta[theta][r] = zeta.imag

"""風圧中心"""
def wind_center():

    """Cl・Cdパラメータ"""
    xp = 0.0
    yp = 0.0
    sxcp = 0.0

    """表面で積分"""
    for theta in range(max_theta - 1):
        dxw = xi[theta + 1][0] - xi[theta][0]
        dyw = eta[theta + 1][0] - eta[theta][0]

        """極座標の方向性により値をマイナスに"""
        dnx =-dyw
        dny = dxw

        """微小部の圧力を端点の平均値とする"""
        cpm = (Cp0[theta + 1][0] + Cp0[theta][0])/2.0

        """微小部の力"""
        fx = cpm*dnx
        fy = cpm*dny

        """微小部の力を足し合わせる"""
        xp = xp + fx #機体軸x方向
        yp = yp + fy #機体軸y方向

        """モーメントの計算のため足し合わせる"""
        sxcp = sxcp + (xi[theta][0]*fy-eta[theta][0]*fx)

    """風圧中心"""
    xcp = sxcp/yp

    """無次元化"""
    cxp = xp / l
    cyp = yp / l

    """機体軸を安定軸方向に"""
    cdp = cxp*math.cos(alfa) + cyp*math.sin(alfa)
    clp = cyp*math.cos(alfa) - cxp*math.sin(alfa)

    """各種係数"""
    print("揚力係数Clは",clp)
    print("抵抗係数Cdは",cdp)
    print("風圧中心xcpは",xcp)

"""等圧線"""
def Cp_plot():

    """範囲指定"""
    plt.xlim(-3.0,3.0)
    plt.ylim(-3.0,3.0)

    """翼を描く(翼はξ、ηの配列のr=0の部分)"""
    plt.plot(xi[:,0],eta[:,0])

    """等圧線を描く"""
    plt.contour(xi[:,:],eta[:,:],Cp0[:,:],locator = plt.MultipleLocator(0.1))
    plt.winter()

    plt.title('Isobar')
    plt.xlabel("ξ")
    plt.ylabel("η")
    cbar = plt.colorbar()
    cbar.set_label("Cp")

    plt.show()

"""表面圧力"""
def Cp_surface():

    """範囲指定"""
    plt.xlim(-1.5,1.5)

    """翼を描く(翼はξ、ηの配列のr=0の部分)"""
    plt.plot(xi[:,0],eta[:,0])

    """表面のCpの値を取る(表面はr=0)"""
    plt.plot(xi[:,0],Cp0[:,0])

    plt.title('Surface Pressure')
    plt.xlabel("ξ")
    plt.ylabel("η")
    plt.show()

if __name__ == "__main__":
    polar()
    Cp_plot()
    Cp_surface()
    wind_center()
