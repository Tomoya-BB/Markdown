import numpy as np
import matplotlib.pyplot as plt

def plot_chord_vs_arc(R=1.0, angle_deg=20):
    """
    2章の図を作りたい。
    Pythonのmatplotlibを用いて、弧(arc)と弦(chord)の差を可視化します。
    
    Parameters
    ----------
    R : float
        球や円の半径 (単位任意。地球半径の縮尺モデルとみなせる)
    angle_deg : float
        P(原点からの距離0の地点)とQ(任意地点)との中心角(度)。
        例: 20度にすると比較的見やすい。
    
    Returns
    -------
    None
        matplotlibの描画を表示
    """
    # 角度をラジアンに変換
    angle_rad = np.radians(angle_deg)
    
    # 円全体を描くための角度(0～2πまで)
    theta_full = np.linspace(0, 2*np.pi, 200)

    # 円(球の断面)のx,y座標 (中心Oは(0,0)に置く)
    x_circle = R * np.cos(theta_full)
    y_circle = R * np.sin(theta_full)

    # 原点Pを例えば、右端(角度0)に置く (便宜上)
    #  -> P = (R,0)
    # Qをangle_deg(ラジアン)だけ回転させた位置 -> Q = (R cos angle, R sin angle)
    # arcはtheta=0からtheta=angle_radまでの円弧
    
    # 円弧描画用(0からangle_radまで)
    theta_arc = np.linspace(0, angle_rad, 100)
    x_arc = R * np.cos(theta_arc)
    y_arc = R * np.sin(theta_arc)
    
    # chord(弦)の始点と終点
    #  P = (R,0) 
    #  Q = (R cos(angle_rad), R sin(angle_rad))
    Px, Py = R, 0
    Qx = R * np.cos(angle_rad)
    Qy = R * np.sin(angle_rad)
    
    # chordを直線として (Px,Py)から(Qx,Qy) への点を生成
    # t in [0,1] で内分
    t_line = np.linspace(0,1,50)
    chord_x = Px + (Qx - Px)*t_line
    chord_y = Py + (Qy - Py)*t_line
    
    # figure, axes
    fig, ax = plt.subplots(figsize=(6,6))
    ax.set_aspect('equal', 'box')  # 1:1のアスペクト比
    
    # 円全体
    ax.plot(x_circle, y_circle, color='lightgray', label='Circle (R=%.2f)' % R)
    
    # 円弧(Arc)
    ax.plot(x_arc, y_arc, color='red', lw=2, label='Arc(0 -> %.1f°)' % angle_deg)
    
    # 弦(Chord)
    ax.plot(chord_x, chord_y, color='blue', lw=2, label='Chord')
    
    # P, Q, O点を描く
    # O(原点)
    ax.plot(0,0, 'ko', label='O(中心)')
    # P
    ax.plot(Px, Py, 'ro', label='P(角度0°)')
    # Q
    ax.plot(Qx, Qy, 'bo', label='Q(角度%.1f°)' % angle_deg)
    
    # 枠やラベル
    ax.set_xlim(-1.2*R, 1.2*R)
    ax.set_ylim(-1.2*R, 1.2*R)
    ax.set_xlabel("X", fontsize=12)
    ax.set_ylabel("Y", fontsize=12)
    ax.set_title("Arc vs Chord (Angle=%.1f°, R=%.1f)" % (angle_deg, R))
    ax.grid(True, alpha=0.3)
    ax.legend(loc='upper right')
    
    plt.show()

# 使い方の例:
if __name__ == "__main__":
    # 半径R=1.0、角度20°の例をプロット
    plot_chord_vs_arc(R=1.0, angle_deg=20)