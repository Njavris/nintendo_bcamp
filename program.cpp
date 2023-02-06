#include <stdio.h>
#include <string.h>
#include <stdint.h>

#define SZ	(sizeof(uint32_t) * 8)
#define MAX_SZ	(256 / SZ)

void encode(int sz, uint32_t *in, uint32_t *out) {
	memset(out, 0, (sz / 16) * sizeof(uint32_t));
	for (int i = 0; i < sz; i++)
		for (int j = 0; j < sz; j++)
			out[(i + j) / 32] ^= ( (in[i / 32] >> (i % 32)) &
				(in[j / 32 + sz / 32] >> (j % 32)) & 1 ) << ((i + j) % 32);
}

struct test {
	int sz;
	uint32_t pln[16];
	uint32_t enc[16];
} tests[] = {
	{ 32, { 0xb0c152f9, 0xebf2831f }, { 0x46508fb7, 0x6677e201 }, },
	{ 32, { 0x00000001, 0x000073af }, { 0x000073af, 0x00000000 }, },
	{ 32, { 0x00000001, 0x738377c1 }, { 0x738377c1, 0x00000000 }, },
	{ 32, { 0xb0c152f9, 0xebf2831f }, { 0x46508fb7, 0x6677e201 }, },
	{ 64, { 0x0cf5c2bf, 0x9aba68ef, 0xc18fb79b, 0xde70eef7 }, { 0xf3268b49, 0x661859eb, 0x0b324559, 0x65ee6bda }, },
	{ 128, { 0xa30d28bd, 0xbda19675, 0x3f95d074, 0xb6f69434, 0xc58f4047, 0xd73fe36a, 0x24be2846, 0xe2ebe432 },
		{ 0xa91db473, 0xfcea8db4, 0xf3bb434a, 0x8dba2f16, 0x51abc87e, 0x92c44759, 0x5c1a16d3, 0x6111c6f4 }, },
};

void dump(int sz, uint32_t *val, uint32_t *exp) {
	printf("%10s","Expected: ");
	for (int i = 0; i < sz / 16; i++)
		printf("%08x ", exp[i]);
	printf("\n%10s","Got: ");
	for (int i = 0; i < sz / 16; i++)
		printf("%08x ", val[i]);
	printf("\n");
}

void run_tests(void) {
	for (int i = 0; i < sizeof(tests) / sizeof(struct test); i++) {
		uint32_t tmp[16];
		struct test *t = &tests[i];
		encode(t->sz, t->pln, tmp);

		if (memcmp(tmp, t->enc, (t->sz / 16) * sizeof(uint32_t)))
			dump(t->sz, tmp, t->enc);
	}
}

void pr_bits(int sz, const uint32_t *poly) {
	for (int i = sz / 32 - 1; i >= 0; i--)
		for (int j = sizeof(uint32_t) * 8 - 1; j >= 0; j--)
			printf("%d", (poly[i] >> j) & 1);
	printf("\n");
}

void pr_poly(int sz, const uint32_t *poly) {
	int empty = 1;
	for (int i = sz / SZ - 1; i >= 0; i--)
		for (int j = sizeof(uint32_t) * 8 - 1; j >= 0; j--)
			if ((poly[i] >> j) & 1) {
				printf("x^%d ", j + 32 * i);
				if ((j + 32 * i) != 0) printf("+ ");
				empty = 0;
			}
	if (empty) printf("0\n");
	else printf("\n");
}

static inline void cpy(int sz, uint32_t *d, const uint32_t *s) {
	for (int i = 0; i < sz / SZ; i++) d[i] = s[i];
}

static inline int pol_ord(int sz, const uint32_t *p) {
	for (int i = sz / SZ - 1; i >= 0; i--)
		for (int j = 31; j >= 0; j--)
			if ((p[i] >> j) & 1)
				return j + 32 * i;
	return -1;
}

static inline void shr(int sz, uint32_t *poly, int bits) {
	uint32_t r = 0;
	uint32_t m = (1 << (bits + 1)) - 1;
	for (int i = sz / SZ - 1; i >= 0; i--) {
		uint32_t t = poly[i] & m;
		poly[i] >>= bits;
		poly[i] |= r << (32 - bits);
		r = t;
	}
}

static inline void shl(int sz, uint32_t *poly, int bits) {
	uint32_t r = 0;
	uint32_t m = ~(1 << bits + 1);
	for (int i = 0; i < sz / SZ; i++) {
		uint32_t t = poly[i] & m;
		poly[i] <<= bits;
		poly[i] |= r >> (32 - bits);
		r = t;
	}
}

static inline void bit(int sz, uint32_t *poly, int bit) {
	poly[bit / SZ] |= 1 << (bit % SZ);
}

static inline int opGF2(int sz, const uint32_t *p1, const uint32_t *p2, uint32_t *res) {
	uint32_t t1[MAX_SZ], t2[MAX_SZ];
	int ret = 0;
	int o1 = pol_ord(sz, p1);
	int o2 = pol_ord(sz, p2);
	cpy(sz, t1, p1);
	cpy(sz, t2, p2);

	if (o2 > o1)
		shl(sz, t1, o2 - o1);
	else if (o1 > o2)
		shl(sz, t2, o1 - o2);

	memset(res, 0, (sz / 32) * sizeof(uint32_t));
	for (int i = 0; i < sz / SZ; i++)
		res[i] = t1[i] ^ t2[i];
	return o1 - o2;
}

void egcd(int sz, const uint32_t *f, uint32_t *res) {
	uint32_t tmp[3][MAX_SZ];
	uint32_t *t1 = tmp[0], *t2 = tmp[1], *t3 = tmp[2];
	memset(res, 0, (sz / SZ) * sizeof(uint32_t));

	/* derivative */
	cpy(sz, t1, f);
	for (int i = sz / SZ - 1; i >= 0; i--)
		t1[i] = (t1[i] & (~0x55555555)) >> 1;

	opGF2(sz, f, t1, t2);
#ifdef DBG
	printf("ECD:\n");
	pr_poly(sz, f);
	pr_poly(sz, t1);
	printf("-----------------------------------\n");
	pr_poly(sz, t2);
#endif
	while (!(pol_ord(sz, t2) < 0)) {
		uint32_t *t = t1;
		t1 = t2;
		t2 = t3;
		t3 = t;
		opGF2(sz, t1, t3, t2);
#ifdef DBG
		pr_poly(sz, t3);
		printf("-----------------------------------\n");
		pr_poly(sz, t2);
#endif
	}
	cpy(sz, res, t1);
}

void div(int sz, const uint32_t *divnd, const uint32_t *divor, uint32_t *quot) {
	uint32_t tmp[2][MAX_SZ];
	uint32_t *t1 = tmp[0], *t2 = tmp[1];
	memset(quot, 0, (sz / SZ) * sizeof(uint32_t));
	cpy(sz, t1, divnd);
#ifdef DBG
	printf("Division:\n");
#endif
	do {
		int r;
		uint32_t *t = t2;
		t2 = t1;
		t1 = t;
		r = opGF2(sz, t2, divor, t1);
		quot[r / SZ] |= 1 << (r % SZ);
#ifdef DBG
		pr_poly(sz, divor);
		printf("----------------------------------- %d\n", r);
		pr_poly(sz, t1);
#endif
	} while (!(pol_ord(sz, t1) < 0));
};

int main() {
	run_tests();
//	uint32_t output[MAX_SZ];
//	uint32_t input[] = { 0xf, 0x1 };
//	encode(32, input, output);
//	pr_poly(32, &input[0]);
//	printf("*\n");
//	pr_poly(32, &input[1]);
//	pr_poly(64, output);
//	printf("%08x %08x\n", output[0], output[1]);

	uint32_t gcd[2], f[2] = { 0x00027fb3 /* 0x17 */, 0 }, q[2];
	pr_poly(64, f);
	egcd(64, f, gcd);
	pr_poly(64, gcd);

	div(64, f, gcd, q);
	pr_poly(64, q);


	return 0;
}
