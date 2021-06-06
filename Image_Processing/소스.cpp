#define _CRT_SECURE_NO_WARNINGS
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
//��Ʈ�� �׷��̽����Ͽ�(��� ��Ʈ�� ����)
//��Ʈ�� ������ ���� �Է�
#define X_SIZE 340
#define Y_SIZE 224
#define PI 3.141592
//�Ǽ� ��� ���� ����ü ����
typedef struct reim {
	double re[Y_SIZE][X_SIZE];
	double im[Y_SIZE][X_SIZE];
}ri;
//�Լ� �켱 ����
int dft2(double arr[][X_SIZE], double im[][X_SIZE], int sel); //Ǫ���� ��ȯ
double* rearr(double arr[][X_SIZE]); //��迭
double* filter(double re[][X_SIZE], int sel); //����
double* moire(double arr[][X_SIZE], double theta); //���Ʒ�
int outf(double arr[][X_SIZE], const char* name, int fsize, int freq); //���
static ri stat; //���� ����ü
char* A;
int pad = 0;
//�Լ� �� �迭 ����� �޸� ������ �ܺο� �������� ����
//
//main�Լ��� �Ǽ�, ����θ� �ٷ�� �迭
double main_re[Y_SIZE][X_SIZE] = { 0, };
double main_im[Y_SIZE][X_SIZE] = { 0, };
//Ǫ���� ��ȯ �Լ� �迭�� ����ü
ri dft_arr;
//��迭 �Լ��� ���Ͽ� static �迭
static double re_arr[Y_SIZE][X_SIZE];
//���� ���Ͽ� �迭
double filt_arr[Y_SIZE][X_SIZE];
//���Ʒ� ���Ͽ� �迭
double moire_arr[Y_SIZE][X_SIZE];
//�̹��� ��� �Լ��� ���� ���Ͽ� �迭
double out_arr[Y_SIZE][X_SIZE];
//
void main() {
	FILE* fpt = fopen("test.bmp", "rb"); //��Ʈ�� ���� ���� ���̳ʸ���
	fseek(fpt, 0, SEEK_END); //��Ʈ�� ������ �� ã��
	int fsize = ftell(fpt); //������ �� �ڸ� = ���ϻ�����
	rewind(fpt); //���� ó������ ��
	A = (char*)malloc(fsize);
	fread(A, sizeof(char), fsize, fpt);
	int xc, yc;
	if (A[18] < 0) xc = 256 + A[18];
	else xc = A[18];
	xc += A[19] * 256;
	printf("width : %d, ", xc); //��Ʈ�� 18, 19�� = ���ϳʺ�
	if (A[22] < 0) yc = 256 + A[22];
	else yc = A[22];
	yc += A[23] * 256;
	printf("height : %d\n", yc); //��Ʈ�� 22, 23�� = ���ϳ���
	if (xc != X_SIZE || yc != Y_SIZE) {
		printf("error : ���� ������ ���Է�, ���α׷� ����\n");
		exit(1);
	}
	fclose(fpt); //���ϴݱ�
	int i, j; //�ݺ�����
	double conv;
	int x = 0, y = 0;
	//54��°���� ���� ����, (R,G,B) �׾ƿø��� �����̰� 0���� �� �����ϹǷ�
	//�׷��̽������� 0.299R + 0.587G + 0.114B = Y
	if (X_SIZE % 4 != 0) pad = 4 - X_SIZE * 3 % 4; //m�� 4�ǹ�� bmp �е�
	for (i = 54; i < fsize; i++) {
		if ((i - 54) % (3 * X_SIZE + pad) < 3 * X_SIZE) {
			if (((i - 54) % (3 * X_SIZE + pad)) % 3 == 0) conv = 0.114;
			else if (((i - 54) % (3 * X_SIZE + pad)) % 3 == 1) conv = 0.587;
			else if (((i - 54) % (3 * X_SIZE + pad)) % 3 == 2) conv = 0.299;
			if (A[i] < 0) main_re[Y_SIZE - 1 - (i - 54) / (3 * X_SIZE + pad)][((i - 54) % (3 * X_SIZE + pad)) / 3] += conv * (256 + A[i]);
			else main_re[Y_SIZE - 1 - (i - 54) / (3 * X_SIZE + pad)][((i - 54) % (3 * X_SIZE + pad)) / 3] += conv * A[i];
		}
	}
	//�׷��̽����� ���Ϸ� ���
	outf(main_re, "1_grayscale.bmp", fsize, 0);
	//���� �߰�
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
	//Ǫ���� ��ȯ �Լ� ����
	dft2(main_re, main_im, 0);
	//�迭 ��迭 �Լ��� �ѱ�� �ٽ� �ޱ�
	temp = rearr(stat.re);
	for (int j = 0; j < Y_SIZE; j++) {
		for (int i = 0; i < X_SIZE; i++) stat.re[j][i] = *(temp + j * X_SIZE + i);
	}
	temp = rearr(stat.im);
	for (int j = 0; j < Y_SIZE; j++) {
		for (int i = 0; i < X_SIZE; i++) stat.im[j][i] = *(temp + j * X_SIZE + i);
	}
	//magnitude ���
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
	//����
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
	//��迭
	temp = rearr(stat.re);
	for (int j = 0; j < Y_SIZE; j++) {
		for (int i = 0; i < X_SIZE; i++) stat.re[j][i] = *(temp + j * X_SIZE + i);
	}
	temp = rearr(stat.im);
	for (int j = 0; j < Y_SIZE; j++) {
		for (int i = 0; i < X_SIZE; i++) stat.im[j][i] = *(temp + j * X_SIZE + i);
	}
	//�� Ǫ����
	dft2(stat.re, stat.im, 1);
	//��� �̹��� ���
	outf(stat.re, "4_final.bmp", fsize, 0);
}
//2���� Ǫ���� ��ȯ
int dft2(double re[][X_SIZE], double im[][X_SIZE], int sel) {
	int u, v, x, y;
	//������ �Ʒ��� (��) ���� ���ٺ��� ���پ� Ǫ����
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
	//���ʿ��� ���������� (��) �� ���ٺ��� ���پ� Ǫ����
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
//��迭 �Լ�
double* rearr(double arr[][X_SIZE]) {
	int n = Y_SIZE / 2, m = X_SIZE / 2, i, j;
	//2��и�
	for (j = 0; j < Y_SIZE - n; j++) {
		for (i = 0; i < X_SIZE - m; i++) {
			re_arr[j][i] = arr[Y_SIZE - 1 - j - n][X_SIZE - 1 - i - m];
		}
	}
	//4��и�
	for (j = Y_SIZE - n; j < Y_SIZE; j++) {
		for (i = X_SIZE - m; i < X_SIZE; i++) {
			re_arr[j][i] = arr[Y_SIZE - 1 - (j - Y_SIZE + n)][X_SIZE - 1 - (i - X_SIZE + m)];
		}
	}
	//1��и�
	for (j = 0; j < Y_SIZE - n; j++) {
		for (i = X_SIZE - m; i < X_SIZE; i++) {
			re_arr[j][i] = arr[Y_SIZE - 1 - j - n][X_SIZE - 1 - (i - X_SIZE + m)];
		}
	}
	//3��и�
	for (j = Y_SIZE - n; j < Y_SIZE; j++) {
		for (i = 0; i < X_SIZE - m; i++) {
			re_arr[j][i] = arr[Y_SIZE - 1 - (j - Y_SIZE + n)][X_SIZE - 1 - i - m];
		}
	}
	return (double*)re_arr;
}
//����
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
//���Ʒ� ���� �߰�
double* moire(double arr[][X_SIZE], double theta) {
	int i, j, p, q, r, thickness, interval, sel; //0�Ͻ� �� ���Ʒ�, 1�Ͻ� �� ���Ʒ� ����
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
//���� ��¿�
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