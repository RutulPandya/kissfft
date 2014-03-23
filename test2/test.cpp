#include "kiss_fft.h"
#include <cstdio>
#include <cstdlib>
#include <vector>
#include <string>

static std::string toString(double f)
{
	char buf[1024];
	sprintf(buf, "%.6f", f);
	while (strlen(buf) > 0 && buf[strlen(buf) - 1] == '0')
		buf[strlen(buf) - 1] = 0;
	if (strlen(buf) > 0 && buf[strlen(buf) - 1] == '.')
		buf[strlen(buf) - 1] = 0;
	return buf;
}

short FFT(short int dir,long m,double *x,double *y)
{
   long n,i,i1,j,k,i2,l,l1,l2;
   double c1,c2,tx,ty,t1,t2,u1,u2,z;

   /* Calculate the number of points */
   n = 1;
   for (i=0;i<m;i++) 
      n *= 2;

   /* Do the bit reversal */
   i2 = n >> 1;
   j = 0;
   for (i=0;i<n-1;i++) {
      if (i < j) {
         tx = x[i];
         ty = y[i];
         x[i] = x[j];
         y[i] = y[j];
         x[j] = tx;
         y[j] = ty;
      }
      k = i2;
      while (k <= j) {
         j -= k;
         k >>= 1;
      }
      j += k;
   }

   /* Compute the FFT */
   c1 = -1.0; 
   c2 = 0.0;
   l2 = 1;
   for (l=0;l<m;l++) {
      l1 = l2;
      l2 <<= 1;
      u1 = 1.0; 
      u2 = 0.0;
      for (j=0;j<l1;j++) {
         for (i=j;i<n;i+=l2) {
            i1 = i + l1;
            t1 = u1 * x[i1] - u2 * y[i1];
            t2 = u1 * y[i1] + u2 * x[i1];
            x[i1] = x[i] - t1; 
            y[i1] = y[i] - t2;
            x[i] += t1;
            y[i] += t2;
         }
         z =  u1 * c1 - u2 * c2;
         u2 = u1 * c2 + u2 * c1;
         u1 = z;
      }
      c2 = sqrt((1.0 - c1) / 2.0);
      if (dir == 1) 
         c2 = -c2;
      c1 = sqrt((1.0 + c1) / 2.0);
   }

   /* Scaling for forward transform */
   if (dir == 1) {
      for (i=0;i<n;i++) {
         x[i] /= n;
         y[i] /= n;
      }
   }
   
   return 1;
}

int main()
{
	std::vector<kiss_fft_cpx> in, out, correct;

	FILE * f = fopen("fft.in", "r");
	double r, i;
	while (!feof(f) && fscanf(f, "%lf", &r) > 0)
	{
		kiss_fft_cpx q;
		q.r = r;
		q.i = 0;
		in.push_back(q);
	}
	fclose(f);

	f = fopen("fft.out", "r");
	while (!feof(f))
	{
		i = 0;
		if (fscanf(f, "%lf", &r) <= 0)
			break;
		int ch = getc(f);
		if (ch != EOF)
		{
			ungetc(ch, f);
			if (ch == '+' || ch == '-')
			{
				fscanf(f, "%lf", &i);
				ch = getc(f);
				if (ch != 'i')
					abort();
			}
		}
		kiss_fft_cpx v;
		v.r = r;
		v.i = i;
		correct.push_back(v);
	}
	fclose(f);

	while (in.size() < 1024)
	{
		kiss_fft_cpx v;
		v.r = 0;
		v.i = 0;
		in.push_back(v);
	}

	out.resize(in.size());

//	memcpy(out.data(), in.data(), in.size() * sizeof(kiss_fft_cpx));
//	FFT(1, out.size(), &out[0].r, &out[0].i);
	kiss_fft_cfg cfg = kiss_fft_alloc((int)in.size(), true, 0, 0);
	kiss_fft(cfg, in.data(), out.data());
	free(cfg);

	for (size_t i = 0; i < in.size(); i++)
	{
//		if (it.i == 0)
//			printf("%s ", toString(it.r).c_str());
//		else
//			printf("%s%s%si ", toString(it.r).c_str(), (it.i > 0 ? "+" : ""), toString(it.i).c_str());

		printf("%10.3f %+-10.3f | %10.3f %+-10.3f\n",
			out[i].r, out[i].i,
			correct[i].r, correct[i].i
		);
	}
}
