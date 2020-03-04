#include <stdio.h>
#include <time.h>
#include <stdint.h>
#include <math.h>
#include <iostream>
#include <vector>

#define pierre_dussart_min 355991

using namespace std;

void print_stats(string algorithm, uint64_t counter, clock_t time){
	algorithm += ':';
	cout << algorithm << string(16 - algorithm.length(), ' ')<< counter;
	printf(" primes in %lf\n", ((double)time)/CLOCKS_PER_SEC);
}

uint64_t pierre_dussart(uint64_t x){
	return((x / log(x)) *( 1 + (1/log(x)) + (2.51/(log(x) * log(x)))));
}

uint64_t prime_estimator(uint64_t i){
	if(i >= pierre_dussart_min)
		return(pierre_dussart(i));
	else
		return (i/2 + 1);
}

void deterministic(uint64_t current, uint64_t range){
	if(current > range)
		return;

	uint64_t counter = 0;
	clock_t time = clock();

	if(range >= 2){
		bool exit_flag = false;
		//preperation
		if(current <= 2 && !exit_flag){
			current = 2;
			counter++;
			if(current < range)
				current++;
			else
				exit_flag = true;
		}
		while(current < 5 && !exit_flag){
			for(uint32_t i = 2; current % i != 0; i++){
				if(i*i > current){
					counter++;
					break;
				}
			}
			if(current < range)
				current++;
			else
				exit_flag = true;
		}
		if(current % 2 == 0 && !exit_flag){
			if(current < range)
				current++;
			else
				exit_flag = true;
		}
		//deterministic algorithm
		while(!exit_flag){
			for(uint32_t i = 3; current % i != 0; i += 2){
				if(i*i > current){
					counter++;
					break;
				}
			}
			if(current +1 < range)
				current += 2;
			else
				exit_flag = true;
		}
	}
	time = clock() - time;
	print_stats("Deterministic", counter, time);
}

void atkin(uint64_t minimum, uint64_t range){
	if(minimum > range)
		return;

	uint64_t current = 0, counter = 0;

	vector<uint64_t> prime_vector(prime_estimator(range), 0);
	prime_vector[0] = 2;


	uint64_t last = 0;
	clock_t time = clock();

	if(range >= 2){
		bool exit_flag = false;
		int i;
		//vector<uint64_t>::iterator it;
		if(current <= 2 && !exit_flag){
			current = 2;
			if(minimum <= 2)
				counter++;
			if(current < range)
				current++;
			else
				exit_flag = true;
		}
		while(current < 5 && !exit_flag){
			for(i = 0;current % prime_vector[i] != 0; i++){
				if(current < prime_vector[i] * prime_vector[i]){
					prime_vector[last+1] = current;
					last ++;
					counter +=  (current >= minimum);
					break;
				}
			}
			if(current < range)
				current++;
			else
				exit_flag = true;
		}
		while(current < minimum && !exit_flag){
			for(i = 0;current % prime_vector[i] != 0; i++){
				if(current < prime_vector[i] * prime_vector[i]){
					prime_vector[last+1] = current;
					last ++;
					break;
				}
			}
			if(current < range)
				current++;
			else
				exit_flag = true;
		}
		if(current % 2 == 0 && !exit_flag){
			if(current < range)
				current++;
			else
				exit_flag = true;
		}
		while(!exit_flag){
			for(i = 0; current % prime_vector[i] != 0; i++){
				if(current < (prime_vector[i] * prime_vector[i])){
					prime_vector[last+1] = current;
					last ++;
					counter ++;
					break;
				}
			}
			if(current +1 < range)
				current += 2;
			else
				exit_flag = true;
		}
	}
	time = clock() - time;
	print_stats("Sieve of Atkin", counter, time);
}

uint64_t modulo(uint64_t p, uint64_t k, uint64_t m){
	if(k == 0)
		return 1;
	else if(k & 1)
		return((p * modulo(p,k-1,m)) % m);
	else
		return(modulo((p*p) % m, k >> 1, m));
}

uint64_t ml_calc(uint8_t a, uint64_t r,uint64_t d,uint64_t current){
	uint64_t x = modulo(a,d,current);
	if(x == 1||x == current - 1)
		return 1;
	for(uint64_t j = 0;j < r-1 && r > 1;j++){
		x = modulo(x, 2, current);
		if (x == current - 1)
			return 1;
	}
	return 0;
}

void miller_rabin(uint64_t current, uint64_t range){
	if(current > range)
		return;

	uint64_t counter = 0;
	uint64_t r, d;

	clock_t time = clock();
	if(range >= 2){
		bool exit_flag = false;
		if(current <= 2 && !exit_flag)
			current = 2;
		if(current == 2){
			counter++;
			if(current < range)
				current++;
			else
				exit_flag = true;
		}
		while(current < 5 && !exit_flag){
			for(uint32_t i = 2; current % i != 0; i++){
				if(i*i > current){
					counter++;
					break;
				}
			}
			if(current < range)
				current++;
			else
				exit_flag = true;
		}
		if(current % 2 == 0 && !exit_flag){
			if(current < range)
				current++;
			else
				exit_flag = true;
		}
		for(; current  <= range && !exit_flag; current += 2){
			if(current == 7 || current == 61){
				counter++;
				continue;
			}
			for(r = 0;((current - 1) >> r) % 2 == 0;r++);
			d = ((current - 1) >> r);
			counter += (ml_calc(2, r, d, current) && ml_calc(7, r, d, current) && ml_calc(61, r, d, current));
		}
	}
	time = clock() - time;
	print_stats("Miller-Rabin", counter, time);
}

int main(int argc, char *argv[]){
	switch(argc){
		case 3:
			cout << "3" << endl;
		case 2:
			cout << "2" << endl;
		case 1:
			cout << "1" << endl;

		default:
			cout << "error" << endl;

	}
	uint64_t min = atoi(argv[1]);
	uint64_t max = atoi(argv[2]);
	printf("\nChecking range [%ld,%ld] for primes\n", min, max);

	cout << string(16, ' ') << prime_estimator(max)  << endl;
	miller_rabin(min, max);
	atkin(min, max);
	deterministic(min, max);

	printf("\n");
	return 0;
}