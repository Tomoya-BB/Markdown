#include <stdio.h>
#include <math.h>

// WGS84定数
#define WGS84_A 6378137.0
#define WGS84_F (1.0/298.257223563)
#define M_PI 3.14159265358979323846

static double deg2rad(double deg) {
    return deg * M_PI / 180.0;
}

// 緯度経度からECEF(X, Y, Z)変換
// lat, lon: degree, h: height [m]
static void latlon_to_ecef(double lat_deg, double lon_deg, double h,
                           double *X, double *Y, double *Z) {
    double e2 = 2.0 * WGS84_F - WGS84_F * WGS84_F;
    double lat = deg2rad(lat_deg);
    double lon = deg2rad(lon_deg);

    double sin_lat = sin(lat);
    double cos_lat = cos(lat);
    double sin_lon = sin(lon);
    double cos_lon = cos(lon);

    double N = WGS84_A / sqrt(1.0 - e2 * sin_lat * sin_lat);

    *X = (N + h) * cos_lat * cos_lon;
    *Y = (N + h) * cos_lat * sin_lon;
    *Z = (N * (1.0 - e2) + h) * sin_lat;
}

// ECEF -> ENU変換
// (X, Y, Z): ECEF座標
// (lat0_deg, lon0_deg, h0): 原点の緯度経度・高さ
// 戻り値: E, N, U (ENU座標)
static void ecef_to_enu(double X, double Y, double Z,
                        double lat0_deg, double lon0_deg, double h0,
                        double *E, double *N, double *U) {
    double X0, Y0, Z0;
    latlon_to_ecef(lat0_deg, lon0_deg, h0, &X0, &Y0, &Z0);

    double dX = X - X0;
    double dY = Y - Y0;
    double dZ = Z - Z0;

    double lat0 = deg2rad(lat0_deg);
    double lon0 = deg2rad(lon0_deg);

    double sin_lat0 = sin(lat0);
    double cos_lat0 = cos(lat0);
    double sin_lon0 = sin(lon0);
    double cos_lon0 = cos(lon0);

    // East軸単位ベクトル
    double e_x = -sin_lon0;
    double e_y =  cos_lon0;
    double e_z = 0.0;

    // North軸単位ベクトル
    double n_x = -sin_lat0 * cos_lon0;
    double n_y = -sin_lat0 * sin_lon0;
    double n_z =  cos_lat0;

    // Up軸単位ベクトル
    double u_x = cos_lat0 * cos_lon0;
    double u_y = cos_lat0 * sin_lon0;
    double u_z = sin_lat0;

    *E = e_x * dX + e_y * dY + e_z * dZ;
    *N = n_x * dX + n_y * dY + n_z * dZ;
    *U = u_x * dX + u_y * dY + u_z * dZ;
}

int main(void) {
    // 原点
    double lat0_deg = 35.0;
    double lon0_deg = 135.0;
    double h0 = 0.0;

    // 対象点(例): 北方向へ約1kmほど離れた地点
    double lat_deg = 35.0001;
    double lon_deg = 135.0;
    double h = 0.0;

    // 対象点をECEFへ
    double X, Y, Z;
    latlon_to_ecef(lat_deg, lon_deg, h, &X, &Y, &Z);

    // ECEF -> ENU
    double E, N, U;
    ecef_to_enu(X, Y, Z, lat0_deg, lon0_deg, h0, &E, &N, &U);

    printf("E = %f [m]\nN = %f [m]\nU = %f [m]\n", E, N, U);

    // 原点からの水平距離
    double distance = sqrt(E*E + N*N);
    printf("水平距離 = %f [m]\n", distance);

    return 0;
}