#include <iostream>
#include <complex>
#include <vector>
#include <string>
#include <time.h>
#include <fstream>
using namespace std;
const static long double PI = 3.1415926535897932384626433832795028841972;

int binReverse(int n,int digit) {
	int result = 0;
	int count = 0;
	vector<int> bin;
	while (n > 0) {
		int r = n % 2;
		bin.push_back(r);
		n = n / 2;
		count++;
	}
	for (int i = 0; i < digit - count; i++) {
		bin.push_back(0);
	}
	int len = bin.size();
	for (size_t i = 0; i < len; i++) {
		result += pow(2, i)*bin[len-i-1];
	}

	return result;
}

void print_fft(vector<complex<double> > input, vector<complex<double> >& output) {
	int N = input.size();
	cout << "DIT-FFT Result:" << endl;
	cout << "Input Series:" << endl;
	for (int i = 0; i < N; i++) {
		cout << input[i] << endl;
	}
	cout << "Output Series:" << endl;
	for (int i = 0; i < N; i++) {
		cout << output[i] << endl;
	}
}

void dft(vector<complex<double> > t, vector<complex<double> >& output) {
	clock_t start, end;
	start = clock();
	//计时开始
	int N = t.size();
	double Theta = -2 * PI / N;
	complex<double> WN(cos(Theta), sin(Theta));

	for (int k = 0; k < N; k++) {
		output[k] = 0;
		for (int n = 0; n < N; n++) {
			output[k] += t[n] * pow(WN, n*k);
		}
	}
	//计时结束
	end = clock();
	double dur = (double)(end - start);
	cout << "DFT Time:" << dur / CLOCKS_PER_SEC << endl;

}

void dit_fft(vector<complex<double> > t, vector<complex<double> >& output) {
	clock_t start, end;
	start = clock();
	//计时开始
	int N = t.size();
	int layer = log2(N);
	vector<int> reverse(N);
	for (int i = 0; i < N; i++) {
		reverse[i] = binReverse(i, layer);
	}
	vector<complex<double> > rotationFactor(N/2);
	double Theta = -2 * PI / N;
	complex<double> WN(cos(Theta), sin(Theta));
	for (int i = 0; i < N / 2; i++) {
		rotationFactor[i] = pow(WN, i);
	}
	for (int i = 0; i < N; i++) {
		output[i] = t[reverse[i]];
	}
	for (int i = 0; i < layer; i++) {
		int DFTsize = pow(2, i + 1);
		int group = N / DFTsize;
		for (int j = 0; j < group; j++) {
			complex<double> temp1, temp2;
			for (int k = 0; k < DFTsize / 2; k++) {
				int n = k + j*DFTsize;
				output[n + DFTsize / 2] = output[n + DFTsize / 2] * rotationFactor[k*group];
				temp1 = output[n] + output[n + DFTsize / 2];
				temp2 = output[n] - output[n + DFTsize / 2];
				output[n] = temp1;
				output[n + DFTsize / 2] = temp2;
			}

		}
	}
	//计时结束
	end = clock();
	double dur = (double)(end - start);
	cout << "DIT-FFT Time:" << dur / CLOCKS_PER_SEC << endl;

}

void dif_fft(vector<complex<double> > t, vector<complex<double> >& output) {
	clock_t start, end;
	start = clock();
	//计时开始
	int N = t.size();
	int layer = log2(N);
	vector<int> reverse(N);
	for (int i = 0; i < N; i++) {
		reverse[i] = binReverse(i, layer);
	}

	vector<complex<double> > rotationFactor(N / 2);
	double Theta = -2 * PI / N;
	complex<double> WN(cos(Theta), sin(Theta));
	for (int i = 0; i < N / 2; i++) {
		rotationFactor[i] = pow(WN, i);
	}
	for (int i = 0; i < N; i++) {
		output[i] = t[i];
	}
	for (int i = 0; i < layer; i++) {
		int group = pow(2, i);
		int DFTsize = N / group;
		for (int j = 0; j < group; j++) {
			complex<double> temp1, temp2;
			for (int k = 0; k < DFTsize / 2; k++) {
				int n = k + j*DFTsize;
				temp1 = output[n] + output[n + DFTsize / 2];
				temp2 = output[n] - output[n + DFTsize / 2];
				output[n] = temp1;
				output[n + DFTsize / 2] = temp2;

				output[n + DFTsize / 2] = output[n + DFTsize / 2] * rotationFactor[k*group];
			}
		}
	}
	vector<complex<double> > temp(N);
	for (int i = 0; i < N; i++) {
		temp[i] = output[reverse[i]];
	}
	for (int i = 0; i < N; i++) {
		output[i] = temp[i];
	}
	//计时结束
	end = clock();
	double dur = (double)(end - start);
	cout << "DIF-FFT Time:" << dur / CLOCKS_PER_SEC << endl;

}

vector<complex<double> > getInput(string filename) {
	ifstream infile(filename);
	vector<complex<double> > xn;
	int fileError = 1;
	try {
		if (!infile) {
			throw fileError;
		}
		string str;
		complex<double> input;
		while (getline(infile, str)) {
			input = stod(str);
			xn.push_back(input);
		}
		return xn;
	}
	catch (int) {
		cout << "FileError:Cannot open file named " <<filename<< endl;
		return xn;
	}
}

int main() {
	string dataset[7] = { "input10.txt","input11.txt","input12.txt","input13.txt","input14.txt","input15.txt","input16.txt" };
	for (int i = 0; i < 7; i++) {
		vector<complex<double> > xn = getInput(dataset[i]);
		int N = xn.size();
		try {
			int digit = log2(N);
			int reverse = binReverse(N - 1, digit);
			if (N - 1 != reverse) {
				throw N;
			}
			cout << "Length:" << N << endl;
			vector<complex<double> > Xk(N);
			vector<complex<double> > Xk1(N);
			vector<complex<double> > Xk2(N);
			dit_fft(xn, Xk);
			dif_fft(xn, Xk1);
			dft(xn, Xk2);
		}

		catch (int) {
			cout << "Input Error:Length of input series must be a power of 2." << endl;
		}
	}

	return 0;
}