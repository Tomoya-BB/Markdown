import math
import tkinter as tk
from tkinter import messagebox
import matplotlib.pyplot as plt

# --------------------------------------------
# 定数および計算用関数
# --------------------------------------------
# WGS84の地球楕円体定数
a = 6378137.0            # 地球長半径 [m]
f = 1.0/298.257223563    # 扁平率
e2 = 2*f - f*f           # 第一離心率^2

def latlon_to_ecef(lat_deg, lon_deg, h=0.0):
    """
    緯度経度(度)と高さh(m)からECEF座標(X,Y,Z)を求める関数。
    この関数により地理座標系から地球中心固定座標系へ変換し、直交座標上でベクトル演算が可能になる。
    """
    lat = math.radians(lat_deg)
    lon = math.radians(lon_deg)
    sin_lat = math.sin(lat)
    cos_lat = math.cos(lat)
    sin_lon = math.sin(lon)
    cos_lon = math.cos(lon)
    
    N = a / math.sqrt(1 - e2*(sin_lat**2))
    X = (N + h) * cos_lat * cos_lon
    Y = (N + h) * cos_lat * sin_lon
    Z = (N*(1 - e2) + h)*sin_lat
    return X, Y, Z

def ecef_to_latlon(X, Y, Z):
    """
    ECEF(X,Y,Z)から緯度経度高さを求める関数。
    Bowring法などを用いてWGS84楕円体上での逆変換を行う。
    """
    # 経度lon
    lon = math.atan2(Y, X)
    
    p = math.sqrt(X*X + Y*Y)
    # 短半径b
    b = a * math.sqrt(1 - e2)
    # 第2離心率^2
    e2p = (a*a - b*b)/(b*b)
    
    # Bowringの方法に基づく計算
    theta = math.atan2(Z * a, p * b)
    sin_theta = math.sin(theta)
    cos_theta = math.cos(theta)
    
    phi = math.atan2(Z + e2p * b * (sin_theta**3),
                     p - e2 * a * (cos_theta**3))
    
    sin_phi = math.sin(phi)
    N = a / math.sqrt(1 - e2 * sin_phi*sin_phi)
    h = p / math.cos(phi) - N
    
    lat_deg = math.degrees(phi)
    lon_deg = math.degrees(lon)
    return lat_deg, lon_deg, h

def compute_projection(A, B, P):
    """
    A,B,P(緯度経度)から、PをA-B直線上に投影した点P'を求める関数。
    手順：
    1. A,B,PをECEF座標へ変換
    2. ベクトルAB,APを算出し、内積計算からパラメータtを求める
    3. P' = A + t*AB をECEFで求める
    4. ECEF→緯度経度逆変換でP'を求める
    """
    X_A, Y_A, Z_A = latlon_to_ecef(A[0], A[1])
    X_B, Y_B, Z_B = latlon_to_ecef(B[0], B[1])
    X_P, Y_P, Z_P = latlon_to_ecef(P[0], P[1])
    
    AB = (X_B - X_A, Y_B - Y_A, Z_B - Z_A)
    AP = (X_P - X_A, Y_P - Y_A, Z_P - Z_A)
    
    APdotAB = AB[0]*AP[0] + AB[1]*AP[1] + AB[2]*AP[2]
    ABdotAB = AB[0]*AB[0] + AB[1]*AB[1] + AB[2]*AB[2]
    
    t = APdotAB / ABdotAB
    
    X_Pp = X_A + t*AB[0]
    Y_Pp = Y_A + t*AB[1]
    Z_Pp = Z_A + t*AB[2]
    
    lat_Pp, lon_Pp, _ = ecef_to_latlon(X_Pp, Y_Pp, Z_Pp)
    return (lat_Pp, lon_Pp)

# --------------------------------------------
# GUI クラス定義
# --------------------------------------------
class App:
    def __init__(self, master):
        master.title("A-Bライン上へのP投影ツール")
        
        # A点入力
        tk.Label(master, text="A 緯度:").grid(row=0, column=0, sticky="e")
        tk.Label(master, text="A 経度:").grid(row=0, column=2, sticky="e")
        self.A_lat = tk.Entry(master)
        self.A_lon = tk.Entry(master)
        self.A_lat.grid(row=0, column=1)
        self.A_lon.grid(row=0, column=3)
        
        # B点入力
        tk.Label(master, text="B 緯度:").grid(row=1, column=0, sticky="e")
        tk.Label(master, text="B 経度:").grid(row=1, column=2, sticky="e")
        self.B_lat = tk.Entry(master)
        self.B_lon = tk.Entry(master)
        self.B_lat.grid(row=1, column=1)
        self.B_lon.grid(row=1, column=3)
        
        # P点入力
        tk.Label(master, text="P 緯度:").grid(row=2, column=0, sticky="e")
        tk.Label(master, text="P 経度:").grid(row=2, column=2, sticky="e")
        self.P_lat = tk.Entry(master)
        self.P_lon = tk.Entry(master)
        self.P_lat.grid(row=2, column=1)
        self.P_lon.grid(row=2, column=3)
        
        # デフォルト値（例）
        self.A_lat.insert(0, "35.394927")
        self.A_lon.insert(0, "136.850982")
        self.B_lat.insert(0, "35.393738")
        self.B_lon.insert(0, "136.888318")
        self.P_lat.insert(0, "35.44375")
        self.P_lon.insert(0, "136.865959")
        
        # 計算＆プロットボタン
        self.compute_button = tk.Button(master, text="計算＆プロット", command=self.compute_and_plot)
        self.compute_button.grid(row=3, column=0, columnspan=4, pady=10)
        
    def compute_and_plot(self):
        # 入力値取得
        try:
            A_lat = float(self.A_lat.get())
            A_lon = float(self.A_lon.get())
            B_lat = float(self.B_lat.get())
            B_lon = float(self.B_lon.get())
            P_lat = float(self.P_lat.get())
            P_lon = float(self.P_lon.get())
        except ValueError:
            messagebox.showerror("入力エラー", "有効な数値を入力してください。")
            return
        
        A = (A_lat, A_lon)
        B = (B_lat, B_lon)
        P = (P_lat, P_lon)
        
        # P'計算
        P_prime = compute_projection(A, B, P)
        
        # 結果表示
        msg = f"P'の計算結果:\n緯度: {P_prime[0]:.6f}\n経度: {P_prime[1]:.6f}"
        messagebox.showinfo("結果", msg)
        
        # プロット: latを縦軸、lonを横軸として単純表示
        plt.figure(figsize=(6,4))
        
        # A,B,P,P'をプロット
        lats = [A[0], B[0], P[0], P_prime[0]]
        lons = [A[1], B[1], P[1], P_prime[1]]
        colors = ["red", "blue", "green", "purple"]
        labels = ["A","B","P","P'"]
        
        for (la, lo, c, lbl) in zip(lats, lons, colors, labels):
            plt.scatter(lo, la, color=c)
            plt.text(lo, la, f" {lbl}", fontsize=9, color=c)
        
        plt.xlabel("経度(度)")
        plt.ylabel("緯度(度)")
        plt.title("A, B, P および P' の位置")
        plt.grid(True)
        plt.show()

# --------------------------------------------
# メインループ開始
# --------------------------------------------
if __name__ == "__main__":
    root = tk.Tk()
    app = App(root)
    root.mainloop()