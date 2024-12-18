#include <stdio.h>
#include <math.h>

// 楕円体定数: GRS80(WGS84とほぼ同等)
#define A 6378137.0
#define F (1.0/298.257222101)
#define E2 (2*F - F*F)
#define DEG2RAD (M_PI/180.0)

//------------------------------------------------------------
// 子午線弧長 M(φ)を計算する関数
// φ(ラジアン)までの子午線弧長を求める。
// 数式はガウス・クリューゲル投影で標準的な展開式。
//------------------------------------------------------------
static double meridian_arc(double phi) {
    double e4 = E2*E2;
    double e6 = e4*E2;

    return A * (
        (1 - E2/4 - 3*e4/64 - 5*e6/256)*phi
        - (3*E2/8 + 3*e4/32 + 45*e6/1024)*sin(2*phi)
        + (15*e4/256 + 45*e6/1024)*sin(4*phi)
        - (35*e6/3072)*sin(6*phi)
    );
}

//------------------------------------------------------------
// ガウス・クリューゲル投影による緯度経度→平面座標変換
// 入力：lat_deg, lon_deg(度単位の緯度経度)
//       φ0_deg(原点緯度), lam0_deg(原点経度), k0(縮尺), FE, FN(偽東北距)
// 出力：X(m), Y(m)
//
// 数式:
// X = FE + k0*N(φ)[A + (1-T+C)A³/6 + ... ]
// Y = FN + k0[M(φ)-M(φ0) + N(φ)*tan(φ)(A²/2 + ...)]
// 詳細な数式はこれまでの解説参照。
//------------------------------------------------------------
static void latlon_to_xy(double lat_deg, double lon_deg,
                         double phi0_deg, double lam0_deg,
                         double k0, double FE, double FN,
                         double *X, double *Y) {
    double phi = lat_deg * DEG2RAD;
    double lam = lon_deg * DEG2RAD;
    double phi0 = phi0_deg * DEG2RAD;
    double lam0 = lam0_deg * DEG2RAD;

    double sinphi = sin(phi);
    double cosphi = cos(phi);
    double N = A / sqrt(1 - E2*sinphi*sinphi);
    double T = (tan(phi)*tan(phi));
    double C = (E2*cosphi*cosphi)/(1-E2);
    double dlam = lam - lam0;
    double A_ = dlam*cosphi;

    double M = meridian_arc(phi);
    double M0 = meridian_arc(phi0);

    double A2 = A_*A_;
    double A3 = A2*A_;
    double A4 = A2*A2;
    double A5 = A4*A_;
    double A6 = A4*A2;

    double x_val = FE + k0*N*(A_ + (1 - T + C)*A3/6.0
        + (5 -18*T + T*T +72*C -58*E2)*A5/120.0);

    double y_val = FN + k0*((M - M0) + N*tan(phi)*(A2/2.0
        + (5 - T +9*C +4*C*C)*A4/24.0
        + (61 -58*T + T*T +600*C -330*E2)*A6/720.0));

    *X = x_val;
    *Y = y_val;
}

int main() {
    double X, Y;

    //------------------------------------------------------------
    // 日本 第7系の例
    // 原点: φ0=36°, λ0=137°10'(=137.1666667°), k0=0.9999, FE=0,FN=0
    // 点: φ=37°N, λ=137.5°E を第7系座標に
    //------------------------------------------------------------
    double lat_jp = 37.0;
    double lon_jp = 137.5;
    latlon_to_xy(lat_jp, lon_jp, 36.0, 137.1666667, 0.9999, 0.0, 0.0, &X, &Y);
    printf("Japan(7th): Lat=%.4f°N Lon=%.4f°E -> X=%.3f m, Y=%.3f m\n",
        lat_jp, lon_jp, X, Y);

    //------------------------------------------------------------
    // デンマーク DKTM各ゾーン
    // DKTM1: φ0=0°, λ0=9°E, k0=1.0, FE=200000m, FN=0
    // DKTM2: φ0=0°, λ0=10°E, k0=1.0, FE=600000m, FN=0
    // DKTM3: φ0=0°, λ0=11°E, k0=1.0, FE=1000000m, FN=0
    //
    // 例: φ=56°N, λ=11°E の点を各DKTMゾーンで計算・比較
    //------------------------------------------------------------
    double lat_dk = 56.0;
    double lon_dk = 11.0;

    // DKTM1
    latlon_to_xy(lat_dk, lon_dk, 0.0, 9.0, 1.0, 200000.0, 0.0, &X, &Y);
    printf("DKTM1: Lat=%.4f°N Lon=%.4f°E -> X=%.3f m, Y=%.3f m\n", lat_dk, lon_dk, X, Y);

    // DKTM2
    latlon_to_xy(lat_dk, lon_dk, 0.0, 10.0, 1.0, 600000.0, 0.0, &X, &Y);
    printf("DKTM2: Lat=%.4f°N Lon=%.4f°E -> X=%.3f m, Y=%.3f m\n", lat_dk, lon_dk, X, Y);

    // DKTM3
    latlon_to_xy(lat_dk, lon_dk, 0.0, 11.0, 1.0, 1000000.0, 0.0, &X, &Y);
    printf("DKTM3: Lat=%.4f°N Lon=%.4f°E -> X=%.3f m, Y=%.3f m\n", lat_dk, lon_dk, X, Y);

    //------------------------------------------------------------
    // 計算例総括:
    // 同じ緯度経度でも、投影パラメータ(φ0, λ0, k0, FE,FN)を変えると
    // X,Y結果が異なり、各地域に最適化した座標系が得られる。
    // 日本とデンマークで使用するパラメータが違うだけで、
    // 数式や計算手順は共通である。
    //------------------------------------------------------------

    return 0;
}