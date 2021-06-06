#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//비트맵 그레이스케일용(흑백 비트맵 파일)
//비트맵 사이즈 정의 입력
#define X_SIZE 340
#define Y_SIZE 224
#define PI 3.141592
//실수 허수 묶음 구조체 선언
typedef struct reim {
	double re[Y_SIZE][X_SIZE];
	double im[Y_SIZE][X_SIZE];
}ri;
//함수 우선 선언
int dft2(double arr[][X_SIZE], double im[][X_SIZE], int sel); //푸리에 변환
double* rearr(double arr[][X_SIZE]); //재배열
double* filter(double re[][X_SIZE], int sel); //필터
double* moire(double arr[][X_SIZE], double theta); //무아레
int outf(double arr[][X_SIZE], const char* name, int fsize, int freq); //출력
static ri stat; //전역 구조체
char* A;
int pad = 0;
//함수 내 배열 선언시 메모리 문제로 외부에 전역변수 선언
//
//main함수의 실수, 허수부를 다루는 배열
double main_re[Y_SIZE][X_SIZE] = { 0, };
double main_im[Y_SIZE][X_SIZE] = { 0, };
//푸리에 변환 함수 배열의 구조체
ri dft_arr;
//재배열 함수의 리턴용 static 배열
static double re_arr[Y_SIZE][X_SIZE];
//필터 리턴용 배열
double filt_arr[Y_SIZE][X_SIZE];
//무아레 리턴용 배열
double moire_arr[Y_SIZE][X_SIZE];
//이미지 출력 함수의 파일 리턴용 배열
double out_arr[Y_SIZE][X_SIZE];
//
void main() {
	FILE* fpt = fopen("test.bmp", "rb"); //비트맵 파일 열기 바이너리로
	fseek(fpt, 0, SEEK_END); //비트맵 파일의 끝 찾기
	int fsize = ftell(fpt); //파일의 끝 자리 = 파일사이즈
	rewind(fpt); //파일 처음으로 감
	A = (char*)malloc(fsize);
	fread(A, sizeof(char), fsize, fpt);
	int xc, yc;
	if (A[18] < 0) xc = 256 + A[18];
	else xc = A[18];
	xc += A[19] * 256;
	printf("width : %d, ", xc); //비트맵 18, 19번 = 파일너비
	if (A[22] < 0) yc = 256 + A[22];
	else yc = A[22];
	yc += A[23] * 256;
	printf("height : %d\n", yc); //비트맵 22, 23번 = 파일높이
	if (xc != X_SIZE || yc != Y_SIZE) {
		printf("error : 정의 사이즈 오입력, 프로그램 종료\n");
		exit(1);
	}
	fclose(fpt); //파일닫기
	int i, j; //반복문용
	double conv;
	int x = 0, y = 0;
	//54번째부터 파일 시작, (R,G,B) 쌓아올리는 형식이고 0으로 줄 구분하므로
	//그레이스케일은 0.299R + 0.587G + 0.114B = Y
	if (X_SIZE % 4 != 0) pad = 4 - X_SIZE * 3 % 4; //m은 4의배수 bmp 패딩
	for (i = 54; i < fsize; i++) {
		if ((i - 54) % (3 * X_SIZE + pad) < 3 * X_SIZE) {
			if (((i - 54) % (3 * X_SIZE + pad)) % 3 == 0) conv = 0.114;
			else if (((i - 54) % (3 * X_SIZE + pad)) % 3 == 1) conv = 0.587;
			else if (((i - 54) % (3 * X_SIZE + pad)) % 3 == 2) conv = 0.299;
			if (A[i] < 0) main_re[Y_SIZE - 1 - (i - 54) / (3 * X_SIZE + pad)][((i - 54) % (3 * X_SIZE + pad)) / 3] += conv * (256 + A[i]);
			else main_re[Y_SIZE - 1 - (i - 54) / (3 * X_SIZE + pad)][((i - 54) % (3 * X_SIZE + pad)) / 3] += conv * A[i];
		}
	}
	//그레이스케일 파일로 출력
	outf(main_re, "1_grayscale.bmp", fsize, 0);
	//잡음 추가
	int sel;
	double* temp;
	printf("\nmoire pattern -> degree n (n <= 90)\nn = -1 none pattren\ninput n : ");
	scanf("%d", &sel);
	if (sel == -1) {}
	else {
		temp = moire(main_re, PI / 180 * sel);
		for (int j = 0; j < Y_SIZE; j++) {
			for (int i = 0; i < X_SIZE; i++) main_re[j][i] = *(temp + j * X_SIZE + i);
		}
		outf(main_re, "1_moire.bmp", fsize, 0);
	}
	//푸리에 변환 함수 실행
	dft2(main_re, main_im, 0);
	//배열 재배열 함수로 넘긴뒤 다시 받기
	temp = rearr(stat.re);
	for (int j = 0; j < Y_SIZE; j++) {
		for (int i = 0; i < X_SIZE; i++) stat.re[j][i] = *(temp + j * X_SIZE + i);
	}
	temp = rearr(stat.im);
	for (int j = 0; j < Y_SIZE; j++) {
		for (int i = 0; i < X_SIZE; i++) stat.im[j][i] = *(temp + j * X_SIZE + i);
	}
	//magnitude 출력
	outf(stat.re, "2_four.bmp", fsize, 1);
	//
	//double mag_avg = 0;
	//double re_avg = 0;
	//double im_avg = 0;
	//for (int j = 0; j < Y_SIZE; j++) {
	//	for (int i = 0; i < X_SIZE; i++) {
	//		mag_avg += out_arr[j][i] / (X_SIZE * Y_SIZE);
	//		re_avg += fabs(stat.re[j][i]) / (X_SIZE * Y_SIZE);
	//		im_avg += fabs(stat.im[j][i]) / (X_SIZE * Y_SIZE);
	//	}
	//}
	//double d0 = X_SIZE / 8;
	//for (int j = 0; j < Y_SIZE; j++) {
	//	for (int i = 0; i < X_SIZE; i++) {
	//		if (out_arr[j][i] - mag_avg * 1.1 > 0 && sqrt(pow(j - Y_SIZE / 2,2) + pow(i - X_SIZE / 2,2)) > d0) {
	//			stat.re[j][i] = stat.re[j][i] / fabs(stat.re[j][i]) * re_avg * 0;// *exp((-1) * (pow(j - Y_SIZE / 2, 2) + pow(i - X_SIZE / 2, 2)) / (2 * d0 * d0));
	//			stat.im[j][i] = stat.im[j][i] / fabs(stat.im[j][i]) * im_avg * 0;// *exp((-1) * (pow(j - Y_SIZE / 2, 2) + pow(i - X_SIZE / 2, 2)) / (2 * d0 * d0));
	//		}
	//	}
	//}
	//필터
	printf("\nfilter select - 0)Ideal, 1)Gaussian, 2)Butterworth, 3)Chebychev 4)Custom\ninput : ");
	scanf("%d", &sel);
	temp = filter(stat.re, sel);
	for (int j = 0; j < Y_SIZE; j++) {
		for (int i = 0; i < X_SIZE; i++) stat.re[j][i] = *(temp + j * X_SIZE + i);
	}
	temp = filter(stat.im, sel);
	for (int j = 0; j < Y_SIZE; j++) {
		for (int i = 0; i < X_SIZE; i++) stat.im[j][i] = *(temp + j * X_SIZE + i);
	}
	outf(stat.re, "3_filtered.bmp", fsize, 1);
	//재배열
	temp = rearr(stat.re);
	for (int j = 0; j < Y_SIZE; j++) {
		for (int i = 0; i < X_SIZE; i++) stat.re[j][i] = *(temp + j * X_SIZE + i);
	}
	temp = rearr(stat.im);
	for (int j = 0; j < Y_SIZE; j++) {
		for (int i = 0; i < X_SIZE; i++) stat.im[j][i] = *(temp + j * X_SIZE + i);
	}
	//역 푸리에
	dft2(stat.re, stat.im, 1);
	//결과 이미지 출력
	outf(stat.re, "4_final.bmp", fsize, 0);
}
//2차원 푸리에 변환
int dft2(double re[][X_SIZE], double im[][X_SIZE], int sel) {
	int u, v, x, y;
	//위에서 아래로 (↓) 왼쪽 끝줄부터 한줄씩 푸리에
	for (u = 0; u < X_SIZE; u++) {
		for (v = 0; v < Y_SIZE; v++) {
			dft_arr.re[v][u] = 0;
			dft_arr.im[v][u] = 0;
			for (y = 0; y < Y_SIZE; y++) {
				dft_arr.re[v][u] = dft_arr.re[v][u] + re[y][u] * cos(2 * PI * v * y / Y_SIZE)
					- im[y][u] * (-1) * sin(2 * PI * v * y / Y_SIZE);
				dft_arr.im[v][u] = dft_arr.im[v][u] + re[y][u] * (-1) * sin(2 * PI * v * y / Y_SIZE)
					- im[y][u] * cos(2 * PI * v * y / Y_SIZE);
			}
		}
	}
	//왼쪽에서 오른쪽으로 (→) 맨 윗줄부터 한줄씩 푸리에
	for (v = 0; v < Y_SIZE; v++) {
		for (u = 0; u < X_SIZE; u++) {
			stat.re[v][u] = 0;
			stat.im[v][u] = 0;
			for (x = 0; x < X_SIZE; x++) {
				stat.re[v][u] = stat.re[v][u] + dft_arr.re[v][x] * cos(2 * PI * u * x / X_SIZE)
					- dft_arr.im[v][x] * (-1) * sin(2 * PI * u * x / X_SIZE);
				stat.im[v][u] = stat.im[v][u] + dft_arr.re[v][x] * (-1) * sin(2 * PI * u * x / X_SIZE)
					- dft_arr.im[v][x] * cos(2 * PI * u * x / X_SIZE);
			}
		}
	}
	if (sel == 1) {
		for (v = 0; v < Y_SIZE; v++) {
			for (u = 0; u < X_SIZE; u++) {
				stat.re[v][u] /= (X_SIZE * Y_SIZE);
			}
		}
	}
	return 0;
}
//재배열 함수
double* rearr(double arr[][X_SIZE]) {
	int n = Y_SIZE / 2, m = X_SIZE / 2, i, j;
	//2사분면
	for (j = 0; j < Y_SIZE - n; j++) {
		for (i = 0; i < X_SIZE - m; i++) {
			re_arr[j][i] = arr[Y_SIZE - 1 - j - n][X_SIZE - 1 - i - m];
		}
	}
	//4사분면
	for (j = Y_SIZE - n; j < Y_SIZE; j++) {
		for (i = X_SIZE - m; i < X_SIZE; i++) {
			re_arr[j][i] = arr[Y_SIZE - 1 - (j - Y_SIZE + n)][X_SIZE - 1 - (i - X_SIZE + m)];
		}
	}
	//1사분면
	for (j = 0; j < Y_SIZE - n; j++) {
		for (i = X_SIZE - m; i < X_SIZE; i++) {
			re_arr[j][i] = arr[Y_SIZE - 1 - j - n][X_SIZE - 1 - (i - X_SIZE + m)];
		}
	}
	//3사분면
	for (j = Y_SIZE - n; j < Y_SIZE; j++) {
		for (i = 0; i < X_SIZE - m; i++) {
			re_arr[j][i] = arr[Y_SIZE - 1 - (j - Y_SIZE + n)][X_SIZE - 1 - i - m];
		}
	}
	return (double*)re_arr;
}
//필터
double* filter(double arr[][X_SIZE], int sel) {
	double x0 = (X_SIZE - 1) / 2;
	double y0 = (Y_SIZE - 1) / 2;
	double e = e = 2.718281;
	int y, x;
	double d;

	double d0 = X_SIZE / 2;
	double n = 2;
	double ep = 1;
	double cn;
	if (sel == 0) {
		for (y = 0; y < Y_SIZE; y++) {
			for (x = 0; x < X_SIZE; x++) {
				d = sqrt((x0 - x) * (x0 - x) + (y0 - y) * (y0 - y));
				if (d > d0) {
					filt_arr[y][x] = 0;
				}
				else {
					filt_arr[y][x] = arr[y][x];
				}
			}
		}
	}
	else if (sel == 1) {
		for (y = 0; y < Y_SIZE; y++) {
			for (x = 0; x < X_SIZE; x++) {
				d = sqrt((x0 - x) * (x0 - x) + (y0 - y) * (y0 - y));
				filt_arr[y][x] = arr[y][x] * exp((-1) * d * d / (2 * d0 * d0));
			}
		}
	}
	else if (sel == 2) {
		for (y = 0; y < Y_SIZE; y++) {
			for (x = 0; x < X_SIZE; x++) {
				d = sqrt((x0 - x) * (x0 - x) + (y0 - y) * (y0 - y));
				filt_arr[y][x] = arr[y][x] * 1 / sqrt(1 + pow((d / d0), 2 * n));
			}
		}
	}
	else if (sel == 3) {
		for (y = 0; y < Y_SIZE; y++) {
			for (x = 0; x < X_SIZE; x++) {
				d = sqrt((x0 - x) * (x0 - x) + (y0 - y) * (y0 - y));
				if (fabs(d / d0) > 1) {
					cn = cosh(n * acosh(d / d0));
				}
				else {
					cn = cos(n * acos(d / d0));
				}
				filt_arr[y][x] = arr[y][x] * 1 / sqrt(1 + ep * ep * cn * cn * (d / d0));
			}
		}
	}
	else if (sel == 4) {
		for (y = 0; y < Y_SIZE; y++) {
			for (x = 0; x < X_SIZE; x++) {
				if (abs(x - 150) < 20) {
					filt_arr[y][x] = arr[y][x];
				}
				else if (abs(x % 50) < 10 && abs(y - 150) < 10) {
					filt_arr[y][x] = 0;
				}
				else {
					filt_arr[y][x] = arr[y][x];
				}
			}
		}
	}
	else {
		for (y = 0; y < Y_SIZE; y++) {
			for (x = 0; x < X_SIZE; x++) {
				filt_arr[y][x] = arr[y][x];
			}
		}
	}
	return (double*)filt_arr;
}
//무아레 잡음 추가
double* moire(double arr[][X_SIZE], double theta) {
	int i, j, p, q, r, thickness, interval, sel; //0일시 선 무아레, 1일시 점 무아레 생성
	printf("\nthickness (integer >= 1) : ");
	scanf("%d", &thickness);
	printf("interval (integer >= thickness) : ");
	scanf("%d", &interval);
	printf("line(0)/dot(1) : ");
	scanf("%d", &sel);
	double tt = tan(theta);
	double a = 0;
	for (j = 0; j < Y_SIZE; j++) {
		for (i = 0; i < X_SIZE; i++) {
			a += arr[j][i] / (X_SIZE * Y_SIZE);
		}
	}
	for (j = 0; j < Y_SIZE; j++) {
		for (i = 0; i < X_SIZE; i++) {
			if (tt < 1) {
				p = (i + 1);
				q = (j + 1) * tt;
				r = (Y_SIZE - j - 1) * tt;
			}
			else {
				p = (i + 1) / tt;
				q = (j + 1);
				r = (Y_SIZE - j - 1);
			}

			if ((q + p) % interval < thickness || (r + p) % interval < thickness) {
				moire_arr[j][i] = 0;
			}
			else {
				moire_arr[j][i] = 255;
			}
		}
	}
	return (double*)moire_arr;
}
//파일 출력용
int outf(double re[][X_SIZE], const char* name, int fsize, int freq) {
	int i, j, l;
	double max = 0;
	if (freq == 1) {
		for (j = 0; j < Y_SIZE; j++) {
			for (i = 0; i < X_SIZE; i++) {
				out_arr[j][i] = log10(sqrt(pow(stat.re[j][i], 2) + pow(stat.im[j][i], 2)));
				if (out_arr[j][i] > max) max = out_arr[j][i];
			}
		}
		for (j = 0; j < Y_SIZE; j++) {
			for (i = 0; i < X_SIZE; i++) {
				if (out_arr[j][i] > 0) {
					out_arr[j][i] = int(out_arr[j][i] * 255 / max);
				}
				if (out_arr[j][i] < 0) {
					out_arr[j][i] = 0;
				}

			}
		}
	}
	else {
		for (j = 0; j < Y_SIZE; j++) {
			for (i = 0; i < X_SIZE; i++) out_arr[j][i] = re[j][i];
		}
	}

	FILE* fpt = fopen(name, "wb");
	fwrite(A, 1, 54, fpt);
	char* B;
	B = (char*)malloc(fsize);
	int k = 0;
	for (j = 0; j < Y_SIZE; j++) {
		for (i = 0; i < X_SIZE; i++) {
			for (l = 0; l < 3; l++) {
				if (out_arr[Y_SIZE - 1 - j][i] < 128) B[k] = out_arr[Y_SIZE - 1 - j][i];
				else B[k] = out_arr[Y_SIZE - 1 - j][i] - 256;
				k++;
			}
		}
		for (l = 0; l < pad; l++) {
			B[k] = 0;
			k++;
		}
	}
	fwrite(B, 1, fsize - 54, fpt);
	fclose(fpt);
	return 0;
}