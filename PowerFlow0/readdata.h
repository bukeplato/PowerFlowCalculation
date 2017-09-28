#include "stdio.h"
#include "math.h"
#include "malloc.h"
#include "string.h"
#include "iostream"
#include "string.h"
using namespace std;
#include "windows.h"
#include "windowsx.h"
#define CRTDBG_MAP_ALLOC
#include "stdlib.h"
#include "conio.h"
#include "crtdbg.h"
/******************************************************全局变量声明*******************************************************/
int node_sum, phd;					//节点数量，平衡节点号
double uph, e;		//平衡点电压，计算精度
int kk;				//作为标识符，标明是否在现有线路基础上叠加，如果叠加就不用开新节点了

/**************************************************函数声明********************************************************/
int readdata();					//数据读取与处理，直角坐标形式，并返回Y阵
//void cY(int);							//转化为幅角形式
char option();
struct Yz *insert1(struct Yz *tp, struct Yz *z);
struct ycb *insert2(struct ycb *, struct ycb *);
/*************************************************结构体定义******************************************************/
struct Ygb		// Y阵（直角坐标）构造
{
	int row, lnxt;			//行列
	double g, b;
	struct Ygb *next;
} *ygb[2000] = { NULL };

struct Ygbd 						//Y阵对角元素（直角坐标）
{
	double g, b;
} ygbd[2000] = { 0 };

struct Yz		// Y阵（幅值相角）构造
{
	int row, lnxt;			//行列
	double y, rad;
	struct Yz *next;
} *y[2000] = { NULL };

struct Yd 						//Y阵对角（幅值相角）元素
{
	double y, rad;
} yd[2000] = { 0 };

struct Ycb
{
	int row, lnxt;
	double zh;
	struct Ycb *next;
} *ycb[4000] = { NULL };

struct pqd
{
	int node;
	double p, q;
	struct pqd *next;
};			//pd点数据
struct pqd *head1;

struct pvd
{
	int node;
	double vi, qmin, qmax;
	struct pvd *next;
};			//pd点数据
struct pvd *head2;

/**************************************************************读取源文件数据***********************************************/
int readdata()
{
	FILE *fp;
	int i, j, m, lastnode = 0;
	char  ch;
	double r, x, bb, e, k0;
	double gij, bij;
	char fsource[20];

	printf("\n		请输入潮流计算的数据文件\n		--");
getf: fflush(stdin);
	gets(fsource);
	while ((fp = fopen(fsource, "rb")) == NULL)
	{
		printf("\t无此文件！  --按[Y]重新输入，其他键退出:");
		fflush(stdin);
		ch = getchar();
		if (ch == 'Y' || ch == 'y')
		{
			printf("\n\t\t潮流计算的数据文件：");
			goto getf;
		}
		else
			printf("\n\t请确认数据文件存在后再重新启动程序！！！");
		Sleep(1000);
		exit(0);
	}

	fscanf(fp, "%d,%d,%lf,%lf", &node_sum, &phd, &uph, &e);		//1、节点数量，平衡节点号，平衡点电压，计算精度
	fscanf(fp, "%d", &m);									//读掉中间的0
	fscanf(fp, "%d", &m);									//读行判断是否有数据
	/************************************************2、除去对角的导纳阵与对角元素********************************************/
	struct Yz *p00, *p11;
	struct Yz *yz_tail = y[1];
	p00 = new Yz;				//直接把他们运算成为幅值相角，三角阵的上半部分
	p11 = new Yz;				//直接把他们运算成为幅值相角，i j颠倒后作为三角阵的下半部分
	while (m != 0)
	{
		fscanf(fp, ",%d,%d,%lf,%lf,%lf", &i, &j, &r, &x, &bb);		// 线路参数
		gij = r / (r*r + x*x);
		bij = -x / (r*r + x*x);


		p00->row = i;
		p00->lnxt = j;
		p00->y = sqrt(gij * gij + bij * bij);
		p00->rad = atan2(bij, gij);
		y[i] = insert1(y[i], p00);
		if (kk == 0)
			p00 = new Yz;				//直接把他们运算成为幅值相角，三角阵的上半部分

		p11->row = j;
		p11->lnxt = i;
		p11->y = sqrt(gij * gij + bij * bij);
		p11->rad = atan2(bij, gij);
		y[j] = insert1(y[j], p11);
		if (kk == 0)
			p11 = new Yz;				//直接把他们运算成为幅值相角，i j颠倒后作为三角阵的下半部分

		ygbd[i].g += gij;
		ygbd[i].b += bij + bb;
		ygbd[j].g += gij;
		ygbd[j].b += bij + bb;
		fscanf(fp, "%d", &m);
		/***************统计y阵链表数组最后的节点号*******************/
		if (i > lastnode)
			lastnode = i;
	}
	delete p00, p11;
	/************************************************3、加上线路中间变压器********************************************/
	fscanf(fp, "%d", &m);
	double x0, b0;
	p00 = new Yz;
	p11 = new Yz;
	while (m != 0)
	{
		fscanf(fp, ",%d,%d,%lf,%lf,%lf", &i, &j, &r, &x, &k0);
		gij = r / (r*r + x*x);
		bij = -x / (r*r + x*x);
		x0 = -gij / k0;
		b0 = -bij / k0;

		p00->row = i;
		p00->lnxt = j;
		p00->y = sqrt(x0 * x0 + b0 * b0);
		p00->rad = atan2(b0, x0);
		y[i] = insert1(y[i], p00);
		if (kk == 0)
			p00 = new Yz;

		p11->row = j;
		p11->lnxt = i;
		p11->y = p00->y;
		p11->rad = p00->rad;
		y[j] = insert1(y[j], p11);
		if (kk == 0)
			p11 = new Yz;

		ygbd[i].g += gij / k0 / k0;
		ygbd[i].b += bij / k0 / k0;
		ygbd[j].g += gij;										//高压侧低压侧有区别
		ygbd[j].b += bij;
		fscanf(fp, "%d", &m);
	}
	delete p00, p11;
	/************************************************4、加上接地支路自导********************************************/
	fscanf(fp, "%d", &m);
	while (m != 0)
	{
		fscanf(fp, ",%d,%lf,%lf", &i, &r, &x);
		ygbd[i].g += r;
		ygbd[i].b += x;
		fscanf(fp, "%d", &m);
	}

	for (int ii = 1; ii <= lastnode; ii++)
	{
		p00 = new Yz;
		p00->row = ii;
		p00->lnxt = ii;
		p00->y = sqrt(ygbd[ii].g * ygbd[ii].g + ygbd[ii].b * ygbd[ii].b);
		p00->rad = atan2(ygbd[ii].b, ygbd[ii].g);
		y[ii] = insert1(y[ii], p00);
	}
	delete p00;

	/************************************************5、读取PQ节点数据********************************************/
	fscanf(fp, "%d", &m);
	struct pqd *p1;
	struct pqd *tail1 = head1;
	while (m != 0)
	{
		fscanf(fp, ",%d,%lf,%lf,%lf,%lf", &i, &r, &x, &bb, &k0);
		p1 = new pqd;
		p1->node = i;
		p1->p = r - bb;
		p1->q = x - k0;

		if (head1 == NULL)
			head1 = p1;
		else
			tail1->next = p1;
		tail1 = p1;
		fscanf(fp, "%d", &m);
	}
	tail1->next = NULL;

	/************************************************6、读取PV节点数据********************************************/
	fscanf(fp, "%d", &m);
	struct pvd *p2;
	struct pvd *tail2 = head2;
	while (m != 0)
	{
		fscanf(fp, ",%d,%lf,%lf,%lf", &i, &r, &x, &bb);
		p2 = new pvd;
		p2->node = i;
		p2->vi = r;
		p2->qmin = x;
		p2->qmax = bb;

		if (head2 == NULL)
			head2 = p2;
		else
			tail2->next = p2;
		tail2 = p2;
		fscanf(fp, "%d", &m);
	}
	tail2->next = NULL;
	fclose(fp);
	/***************统计y阵链表数组最后的节点号*******************/
	return lastnode;
}
