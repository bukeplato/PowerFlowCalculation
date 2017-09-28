#include "stdio.h"
#include "iostream"
using namespace std;

#include "math.h"
#include "string.h"
#include "windows.h"

/******************************************************全局变量声明*******************************************************/
int node_sum, phd;					//节点数量，平衡节点号
double uph, e;						//平衡点电压，计算精度
int kk;								//作为标识符，标明是否在现有线路基础上叠加，如果叠加就不用开新节点了
double v[1500] = { 0 };						//电压幅值
double rad[1500] = { 0 };			//电压相角
double radv[3000] = { 0 };			//奇数为 相角 偶数为 电压幅值
double err = 0;						//误差
int pvyjbz[1500] = { 0 };			//pv点越界标志
double p[1500] = { 0 };			// Pi
double q[1500] = { 0 };			// Qi
double gmin[1500] = { 0 };
double gmax[1500] = { 0 };

/*******************************************************函数声明********************************************************/
double readdata();					//数据读取与处理，直角坐标形式，并返回Y阵
char option();
struct Yz *insert1(struct Yz * , struct Yz * );
struct Ycb *insert2(struct Ycb * , struct Ycb * );
struct Ycb *insert3(struct Ycb * , struct Ycb * );
void ycbmatrix( struct Yz *[] );
void ycbdata( struct Yz *[] );
void luss();
void fixradv();
void deltapq();
void fixphdpvd();
void deletepvd(int);
void result();
void finaldeal();

/*******************************************************结构体定义******************************************************/
struct Ygb		// Y阵（直角坐标）构造
{
	int row, lnxt;			//行列
	double g, b;
	struct Ygb *next;
} *ygb[50] = { NULL };

struct Ygbd 						//Y阵对角元素（直角坐标）
{
	double g, b;
} ygbd[50] = { 0 };

struct Yz		// Y阵（幅值相角）构造
{
	int row, lnxt;			//行列
	double y, rad;
	struct Yz *next;
} *y[50] = { NULL };

struct Yd 						//Y阵对角（幅值相角）元素
{
	double y, rad;
} yd[50] = { 0 };

struct Ycb
{
	int row, lnxt;
	double zh;
	struct Ycb *next;
} *ycb[100] = { NULL };

struct pqd
{
	int node;
	double p, q;
	struct pqd *next;
};			//pd点数据
struct pqd *headpqd;

struct pvd
{
	int node;
	double vi, qmin, qmax;
	struct pvd *next;
};			//pd点数据
struct pvd *headpvd;

/*******************************************************主函数再此处********************************************************/
void main()
{
	char c;						// record input choice
	int ynum = 0;				// 统计y阵链表数组最后的节点号（为了节约空间，只记录了上三角内容，互导对称部分没有记录）

	opt: c = option();			// opt分支指向主界面函数

	if (c == '1')				//键入1时进行潮流运算
		goto opt1;

	else if(c=='2')
	{
	system( "mode con:cols=100 lines=30 & color 07" );

		system("cls");
		printf("         ################################计算结果查询################################\n");
		printf("\n\n");
		result();
		goto opt;
	}

	else if (c == '0')
	{
		char ext;
		printf("\t确定退出？ --按[Y]确认，其他键取消： ");
		fflush(stdin);
		ext = getchar();
		if (ext == 'Y' || ext == 'y')
		{
			printf("\n\t正在退出程序.");
			Sleep(333);
			printf(".");
			Sleep(333);
			printf(".");
			Sleep(333);
			exit(0);
		}
		else
			goto opt;
	}

	else
	{
		printf("\t\t\t【【【指令错误，请重输！】】】\n");
		Sleep(1000);
		goto opt;
	}

	opt1 : e = readdata();

	deltapq();						//不平衡量△P △Q ，存在radv[]，统计err

	int loop = 0;					//循环次数
	while( err > e && loop < 20 )
	{
		err = 0;

		ycbmatrix(y);
		ycbdata(y);
		luss();				// LU 分解得到 △θ △V ，存在radv[]

		fixradv();			// 用上面数据修正后的 θ V ，存在 rad[] v[]
		//fixphdpvd();		//新的 phd、pvd 的 Pi Qi ，存在 p[i] q[i]

		deltapq();			//不平衡量△P △Q ，存在radv[]，统计err
		loop++;
	}


	finaldeal();
/*******************************************文件输出**********************************/
	if(loop >= 20)
	{
		printf("\t\t\t【【【潮流运算不收敛！】】】\n");
		Sleep(1000);
	}

	FILE *fp;
	if((fp=fopen("solution.txt","w")) == NULL)
	{
	printf("文件写入错误！！");
	Sleep(1000);
	exit(0);
	}
	fprintf(fp , "Node Amount:%d\n" , node_sum);
	if(loop >= 20)
	{
		printf("\t\t\t【【【潮流运算不收敛！】】】\n");
		goto end;
	}

	fprintf(fp, "Converged, Iterations:%d\n", loop);
	fprintf(fp, "Node,    Voltage,    Radian,    P,    Q\n");
	fprintf(fp, "-\n");
	for( int i = 1 ; i <= node_sum ; i++ )
		fprintf(fp , "%d,  %f,  %f,  %f,  %f\n" , i , v[i] , rad[i] , p[i] , q[i]);


	printf("\n\t\t\t潮流计算完成!\n");
end:fclose(fp);
	Sleep(500);
	printf("\n\t\t\t文件写入完成!\n");
	printf("\t\t\t按[0]退出程序\n");
	printf("\t\t\t按[b]退回主界面\n");


	while (1)
	{
		c = getchar();
		if (c == '0')
			break;
		if (c == 'b' || c == 'B')
			goto opt;
	}
}



/********************************************* 最后计算phd pvd数值 *************************************/
void finaldeal()
{
	double x = 0, r = 0;
	struct Yz *p1 = NULL;								//平衡点处理
	p1 = y[phd];

	while (p1 != NULL)
	{
		x = x + p1->y * v[p1->lnxt] * cos(rad[phd] - rad[p1->lnxt] - p1->rad);
		r = r + p1->y * v[p1->lnxt] * sin(rad[phd] - rad[p1->lnxt] - p1->rad);
		p1 = p1->next;
	}
	/******************是否要把对角元素单独考虑？**************************/
	x = x + ygbd[phd].g *  v[phd] * cos(rad[phd] - rad[phd] - ygbd[phd].b);
	r = r + ygbd[phd].g *  v[phd] * sin(rad[phd] - rad[phd] - ygbd[phd].b);
	p[phd] = x * v[phd];
	q[phd] = r * v[phd];

	//*****************************PV点处理，p2为 pvd 链表***************************

	struct pvd *p2 = headpvd;
	struct Yz *p3 = NULL;

	while (p2 != NULL)
	{
		r = 0;
		p3 = y[p2->node];		// p2为 pvd 链表, 取其节点号，得到导纳阵表头
		while (p3 != NULL)
		{
			r = r + p3->y * v[p3->lnxt] * sin(rad[p3->row] - rad[p3->lnxt] - p3->rad);
			p3 = p3->next;
		}
		p2 = p2->next;
	}
}

/****************************************************主菜单选项页面*********************************************************/
char option()
{
	char c;
	system("mode con:cols=75 lines=25 &color 70");
	puts("\n");
	puts("      ***************************************************************");
	puts("		          潮流计算课程设计		");
	puts("      ***************************************************************\n\n");
	printf("\t操作指令：\n\n\t\t1.潮流计算      2.计算结果查看      0.退出程序\n\n");
	fflush(stdin);
	printf("\t输入指令：");
	c = getchar();
	return c;
}

/****************************************************结果查询页面*********************************************************/
void result()
{
	int ch;
	char cint = NULL;
	FILE *fp;
	system("mode con:cols=75 lines=25 &color 70");
	puts("\n\n");
	puts("请输入待查节点号：");
	
	cin >> ch;
	
	if (0 < ch && ch <32767)
	{
		
	}
	else
	{
		printf("\t【【【错误的节点号！】】】\n");
		printf("\n\t按[b]返回主界面：");
		while (1)
		{
			cint = NULL;
			cint = getchar();
			if (cint == 'b' || cint == 'B')
				return;
		}
	}

	if((fp = fopen("solution.txt", "rb")) == NULL)
	{
		printf("\n\t运算结果文件丢失！！！");
		Sleep(1000);
		exit(0);
	}


	while (cint != '-')
		fscanf(fp, "%c", &cint);

	printf("\n\n\t节点      电压      相角          P          Q\n\n");

	int node = 0;
	double volt, degree, P, Q;

	while (node != ch )
	{
		fscanf(fp, "%d,  %lf,  %lf,  %lf,  %lf", &node, &volt, &degree, &P, &Q);
	}
	printf("\t%d   || %.5lf || %.5lf || %.5lf || %.5lf", node, volt, degree, P, Q);
	
	printf("\n\n\n\n\t按[0]退出程序：");
	printf("\n\t按[b]返回主界面：");
	while (1)
	{
		cint = getchar();
		if (cint == 'b' || cint == 'B')
			break;
		if (cint == '0' )
			exit(0);
	}

	fclose(fp);
}

/***************************************************读取源文件数据**********************************************************/
double readdata()
{
	FILE *fp;
	int i, j, m;
	char  ch;
	double r, x, bb, ee, k0;
	double gij, bij;
	char fsource[20];

	printf("\n		请输入潮流计算的数据文件\n		--");
	getf: fflush(stdin);
	gets(fsource);
	if ((fp = fopen(fsource, "rb")) == NULL)
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

	fscanf(fp, "%d,%d,%lf,%lf", &node_sum, &phd, &uph, &ee);		//1、节点数量，平衡节点号，平衡点电压，计算精度
	fscanf(fp, "%d", &m);									//读掉中间的0
	fscanf(fp, "%d", &m);									//读行判断是否有数据

	for (int ii = 1; ii <= node_sum; ii++)
		v[ii] = 1;

	/***********************2、除去对角的导纳阵与对角元素**************************/
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
		p00->rad = atan2( - bij , - gij );
		y[i] = insert1(y[i], p00);
		if (kk == 0)
			p00 = new Yz;				//直接把他们运算成为幅值相角，三角阵的上半部分

		p11->row = j;
		p11->lnxt = i;
		p11->y = sqrt(gij * gij + bij * bij);
		p11->rad = atan2(- bij, - gij);
		y[j] = insert1(y[j], p11);
		if (kk == 0)
			p11 = new Yz;				//直接把他们运算成为幅值相角，i j颠倒后作为三角阵的下半部分

		ygbd[i].g += gij;
		ygbd[i].b += bij + bb;
		ygbd[j].g += gij;
		ygbd[j].b += bij + bb;
		fscanf(fp, "%d", &m);
	}
	delete p00, p11;

	/****************************************3、加上线路中间变压器*********************************/
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
		p00->y = sqrt( x0 * x0 + b0 * b0 );
		p00->rad = atan2( b0 , x0 );
		y[i] = insert1( y[i] , p00 );
		if (kk == 0)
			p00 = new Yz;

		p11->row = j;
		p11->lnxt = i;
		p11->y = sqrt(x0 * x0 + b0 * b0);
		p11->rad = atan2(b0, x0);
		y[j] = insert1( y[j] , p11 );
		if (kk == 0)
			p11 = new Yz;

		ygbd[i].g += gij / k0 / k0;
		ygbd[i].b += bij / k0 / k0;
		ygbd[j].g += gij;										//高压侧低压侧有区别
		ygbd[j].b += bij;
		fscanf(fp, "%d", &m);
	}
	delete p00, p11;

	/******************************4、加上接地支路自导***************************/
	fscanf(fp, "%d", &m);
	while (m != 0)
	{
		fscanf(fp, ",%d,%lf,%lf", &i, &r, &x);
		ygbd[i].g += r;
		ygbd[i].b += x;
		fscanf(fp, "%d", &m);
	}

	for (int ii = 1; ii <= node_sum; ii++)
	{
		r = sqrt(ygbd[ii].g * ygbd[ii].g + ygbd[ii].b * ygbd[ii].b);
		ygbd[ii].b = atan2(ygbd[ii].b , ygbd[ii].g);					//变为对角线的相角
		ygbd[ii].g = r;													//变为对角线的幅值
	}

	j = 0;
	for (int ii = 1; ii <= node_sum; ii++)								//判断孤点
	{
		if(ygbd[ii].g == 0)
			j++;
	}
	if(j > 0)
	{
		printf("\n\t\t\t【【【系统中存在孤点！】】】" );
		printf("\n\t\t\t【【【请检查后重启程序！】】】");
		Sleep(3000);
		exit(0);
	}


	/*
	for (int ii = 1; ii <= node_sum; ii++)
	{
		p00 = new Yz;
		p00->row = ii;
		p00->lnxt = ii;
		p00->y = sqrt( ygbd[ii].g * ygbd[ii].g + ygbd[ii].b * ygbd[ii].b );
		p00->rad = atan2( ygbd[ii].b , ygbd[ii].g );
		y[ii] = insert1(y[ii], p00);
	}
	delete p00;
	*/

	/********************************5、读取所有节点功率数据***************************/
	fscanf(fp, "%d", &m);
	/*	struct pqd *p1;
	struct pqd *tailpqd = headpqd;*/
	while (m != 0)
	{
		fscanf(fp, ",%d,%lf,%lf,%lf,%lf", &i, &r, &x, &bb, &k0);

		p[i] = r - bb;						// Pi
		q[i] = x - k0;						// Qi

	/*	p1 = new pqd;
		p1->node = i;
		p1->p = r - bb;
		p1->q = x - k0;

		if (headpqd == NULL)
			headpqd = p1;
		else
			tailpqd->next = p1;
		tailpqd = p1;*/
		fscanf(fp, "%d", &m);
	}
	//tailpqd->next = NULL;

	/********************************6、读取PV节点数据***************************/
	fscanf(fp, "%d", &m);
	struct pvd *p2;
	struct pvd *tailpvd = headpvd;
	while (m != 0)
	{
		fscanf(fp, ",%d,%lf,%lf,%lf", &i, &r, &x, &bb);
		p2 = new pvd;
		p2->node = i;
		p2->vi = r;
		p2->qmin = x;
		p2->qmax = bb;

		gmax[i] = bb;
		gmin[i] = x;
		//pvyjbz[i] = 0;					//pv点越界标志

		v[i] = p2->vi;

		if (headpvd == NULL)
			headpvd = p2;
		else
			tailpvd->next = p2;
		tailpvd = p2;
		fscanf(fp, "%d", &m);
	}
	p[phd] = 0;
	q[phd] = 0;
	radv[2 * phd - 1] = 0;		// phd, △Pi  △Qi = 0
	radv[2 * phd] = 0;
	tailpvd->next = NULL;
	fclose(fp);

	return ee;
}

/***************************************************形成雅可比矩阵稀疏表（竖着排）*******************************************/
void ycbmatrix(struct Yz *tp[ ])
{
	struct Yz *p1;
	struct Ycb *p3, *p4;
	int i, j;

	for( i = 1 ; i <= node_sum ; i++ )		//（只是定位置）
	{
		p1 = tp[i];
		j = 0;								//换行标志
		if( p1 != NULL )
		{
			while(( p1->row == i ) && ( j == 0 ))
			{
				p3 = new Ycb;					//Hij,Nij
				p4 = new Ycb;
				p3->row = p1->row * 2 - 1;
				p4->row = p1->row * 2 - 1;
				p3->lnxt = p1->lnxt * 2 - 1;
				p4->lnxt = p1->lnxt * 2;
				p3->zh = 1;
				p4->zh = 1;
				p3->next = NULL;
				p4->next = NULL;
				ycb[p3->lnxt] = insert2( ycb[p3->lnxt] , p3 );
				ycb[p4->lnxt] = insert2( ycb[p4->lnxt] , p4 );
				if( p1->next == NULL )
					j = 1;
				else
					p1 = p1->next;
			}
		}

		j = 0;									//换行标志
		p1 = tp[i];
		if( p1 != NULL )
		{
			while(( p1->row == i ) && ( j == 0 ))
			{
				p3 = new Ycb;					//Jij,Lij
				p4 = new Ycb;
				p3->row = p4->row = p1->row * 2;
				p3->lnxt = p1->lnxt * 2 - 1;
				p4->lnxt = p1->lnxt * 2;
				p3->zh = 1;
				p4->zh = 1;
				p3->next = NULL;
				p4->next = NULL;
				ycb[p3->lnxt] = insert2( ycb[p3->lnxt] , p3 );
				ycb[p4->lnxt] = insert2( ycb[p4->lnxt] , p4) ;
				if( p1->next == NULL )
					j = 1;
				else
					p1 = p1->next;
			}
		}
		/***************************************对角元素单独处理******************************/
		p3 = new Ycb;							//Hii
		p3->row = 2 * i - 1;
		p3->lnxt = 2 * i - 1;
		p3->zh = 1;
		p3->next = NULL;
		ycb[p3->lnxt] = insert2( ycb[p3->lnxt] , p3 );

		p3 = new Ycb;							//Lii
		p3->row = 2 * i;
		p3->lnxt = 2 * i;
		p3->zh = 1;
		p3->next = NULL;
		ycb[p3->lnxt] = insert2( ycb[p3->lnxt] , p3 );

		p3 = new Ycb;							//Nii
		p3->row = 2 * i - 1;
		p3->lnxt = 2 * i;
		p3->zh = 1;
		p3->next = NULL ;
		ycb[p3->lnxt] = insert2( ycb[p3->lnxt] , p3 );

		p3 = new Ycb;							//Jii
		p3->row = 2 * i;
		p3->lnxt= 2 * i - 1;
		p3->zh = 1;
		p3->next = NULL;
		ycb[p3->lnxt] = insert2( ycb[p3->lnxt] , p3 );
	}
}

/***************************************************形成雅可比矩阵********************************************************/
void ycbdata(struct Yz *tp[ ])
{
	struct Yz *p1;
	struct Ycb *p3;
	int i, j, n, k;
	double x;

	for( i = 1 ; i <= node_sum ; i++ )
	{
		p3 = ycb[ 2 * i - 1 ];				// 对 H J 的列进行处理
		while( p3 != NULL )
		{
			j = p3->row;					// j jcb[]的行号
			n = ( j + 1 ) / 2;						// int n为 节点导纳矩阵 Y的行号，i是列号

			if( fmod( j , 2 ) != 0 )				//奇数，为H
			{
				if( n == i )							//**************** 表示Hii *******************
				{
					p1 = tp[n];
					x = 0;
					k = 0;
					while(( p1->row == n ) && ( k == 0 ))			// 取的行号为n的 y 阵元素值
					{
						x = x + p1->y * v[p1->lnxt] * sin( rad[p1->row] - rad[p1->lnxt] - p1->rad );
						if( p1->next == NULL )
							k = 1;
						else
							p1 = p1->next;
					}
					x = x * v[n];
					p3->zh = x;
				}
				else									//**************** 表示Hij *******************
				{
					p1 = tp[n];
					while( p1->lnxt != i )				// int n为 节点导纳矩阵 Y的行号，i是列号
						p1 = p1->next;
					p3->zh = - v[p1->row] * p1->y * v[p1->lnxt] * sin( rad[p1->row] - rad[p1->lnxt] - p1->rad );
				}
			}

			else								//偶数，为J
			{	if( n == i )						//**************** 表示Jii *******************
				{
					p1 = tp[n];
					x = 0;
					k = 0;
					while(( p1->row == n ) && ( k == 0 ))
					{
						x = x + p1->y * v[p1->lnxt] * cos( rad[p1->row] - rad[p1->lnxt] - p1->rad );
						if( p1->next == NULL )
							k = 1;
						else
							p1 = p1->next;
					}
					x = x * v[n];
					p3->zh = -x;
				}
				else				//**************** 表示Jij  *******************
				{
					p1 = tp[n];
					while( p1->lnxt != i )
						p1 = p1->next;
					p3->zh = v[p1->row] * p1->y * v[p1->lnxt] * cos( rad[p1->row] - rad[p1->lnxt] - p1->rad );
				}
			}
			p3 = p3->next;
		}


		p3 = ycb[ 2 * i ];		//这部分是N或L
		while( p3 != NULL )
		{
			j = p3->row;
			n = ( j + 1 ) / 2;		//n对应的行号，i是列号
			if( fmod( j , 2 ) != 0)		//表示N
			{
				if( n == i )			//表示Nii
				{
					p1 = tp[n];
					x = 0;
					k = 0;
					while(( p1->row == n ) && ( k == 0 ))
					{
						x = x + p1->y * v[p1->lnxt] * cos( rad[p1->row] - rad[p1->lnxt] - p1->rad );
						if( p1->next == NULL )
							k = 1;
						else
							p1 = p1->next;
					}
	/*##################################### 注意 N、L的定义 ################################*/
	/*##################################### 注意 N、L的定义 ################################*/
	/*##################################### 注意 N、L的定义 ################################*/
					//双份的 Vi * Yii * Vi * cos(xita ij)     但是没把附加Vi写入
					p3->zh = - ( x + 2 * v[n] * ygbd[i].g * cos( ygbd[i].b ));
				}
				else				//表示Nij
				{
					p1 = tp[n];
					while( p1->lnxt != i )
						p1 = p1->next;
					p3->zh = - v[p1->row] * p1->y * cos( rad[p1->row] - rad[p1->lnxt] - p1->rad );
				}
			}
			else						//表示L
			{
				if( n == i )			//表示Lii
				{
					p1 = tp[n];
					x = 0;
					k = 0;
					while(( p1->row == n ) && ( k == 0 ))
					{
						x = x + p1->y * v[p1->lnxt] * sin( rad[p1->row] - rad[p1->lnxt] - p1->rad);
						if( p1->next == NULL )
							k = 1;
						else
							p1 = p1->next;
					}
					//双份的 Vi * Yii * Vi * sin(xita ij)
					p3->zh = - ( x - 2 * v[n] * ygbd[i].g * sin( ygbd[i].b ) );
				}
				else				//表示Lij
				{
					p1 = tp[n];
					while( p1->lnxt != i )
						p1 = p1->next;
					p3->zh = - v[p1->row] * p1->y * sin( rad[p1->row] - rad[p1->lnxt] - p1->rad );
				}
			}
			p3 = p3->next;
		}
	}

	struct Ycb  *p4, *p2;
	/****************************************处理平衡节点******************************************************/
	for (i = 1; i <= 2 * node_sum; i++)
	{
		if ((i == 2 * phd) || (i == 2 * phd - 1))	//ycb是按列排的，将整列删除
		{
			while (ycb[i] != NULL)
			{
				p3 = ycb[i]->next;
				delete ycb[i];
				ycb[i] = p3;
			}
		}
		else													//只是删除对应的那一行的元素
		{
			n = 0;
			p3 = ycb[i];
			p4 = ycb[i];
			p2 = ycb[i];
			while (p3 != NULL)
			{
				p4 = p3;
				p3 = p3->next;									// P3 为 P4 的下一个
				if ((p4->row == 2 * phd) || (p4->row == 2 * phd - 1))
				{
					if (n == 0)
					{
						ycb[i] = p3;
						delete p4;
					}
					else
					{
						p2->next = p3;
						delete p4;
					}
				}
				else				//	对应PQ行的元素删光后 n置1
				{
					n = 1;
					p2 = p4;
				}
			}
		}
	}
	/****************************************处理PV节点******************************************************/
	struct pvd *p7 = headpvd;
	while (p7 != NULL)
	{
		for (i = 1; i <= 2 * node_sum; i++)
		{

			if ( i == 2 * p7->node )	//ycb是按列排的，将整列删除	 PV的 Q 对应行列
			{
				while (ycb[i] != NULL)
				{
					p3 = ycb[i]->next;
					delete ycb[i];
					ycb[i] = p3;
				}
			}
			else													//只是删除对应的那一行的元素
			{
				n = 0;
				p3 = ycb[i];
				p4 = ycb[i];
				p2 = ycb[i];
				while (p3 != NULL)
				{
					p4 = p3;
					p3 = p3->next;									// P3 为 P4 的下一个
					if ( p4->row == 2 * p7->node )
					{
						if (n == 0)
						{
							ycb[i] = p3;
							delete p4;
						}
						else
						{
							p2->next = p3;
							delete p4;
						}
					}
					else				//	对应PQ行的元素删光后 n置1
					{
						n = 1;
						p2 = p4;
					}
				}
			}
		}
		p7 = p7->next;
	}

	p3 = new Ycb;				//平衡点要留对角线元素
	p3->row = 2 * phd - 1;
	p3->lnxt = 2* phd - 1;
	p3->zh = 99999;
	p3->next = NULL;
	ycb[p3->lnxt] = p3;

	p3 = new Ycb;				//平衡点要留对角线元素
	p3->row = 2 * phd;
	p3->lnxt  = 2 * phd;
	p3->zh = 99999;
	p3->next = NULL;
	ycb[p3->lnxt] = p3;

	p7 = headpvd;
	while (p7 != NULL)
	{
		p3 = new Ycb;				//PV点 Q 要留对角线元素
		p3->row = 2 * p7->node;
		p3->lnxt = 2 * p7->node;
		p3->zh = 99999;
		p3->next = NULL;
		ycb[p3->lnxt] = p3;
		p7 = p7->next;
	}
}

/********************************************** LU分解求修正方程 ***********************************************************/
void luss()
{
	struct Ycb *l[100] = { NULL }, *u[100] = { NULL };	//分别为雅可比矩阵LU分解稀疏表
	struct Ycb *p1, *p2;
	double ud[100];	//对角线元素
	double t[100];
	double x;
	int i, j, r;

	for(int ii = 1; ii <= 2 * node_sum; ii++)
		radv[ii] = -radv[ii];

	/**********************将雅可比矩阵放到对应的链表里***************************/
	p1 = new Ycb;
	for( i = 1 ; i <= 2 * node_sum ; i++ )
	{
		ud[i] = 0;
		p2 = ycb[i];

		while( p2 != NULL )
		{

			p1->row = p2->row;
			p1->lnxt = p2->lnxt;
			p1->zh = p2->zh;
			p1->next = NULL;

			if( p2->row == p2->lnxt )					//ud 雅克比矩阵的对角值 即为U阵对角
				ud[p2->row] = p2->zh;

			if( p1->row > p1->lnxt )					//l,l矩阵是按行排列的，其行号>列号
			{
				l[p1->row] = insert3( l[p1->row] , p1 );
				if( kk == 0 )							// kk存在于insert函数中，如果有节点间多线路重叠的情况kk为1
					p1 = new Ycb;
			}

			if( p1->row < p1->lnxt )					//u,u矩阵是按列排列的，其列号>行号
			{
				u[p1->lnxt] = insert2( u[p1->lnxt] , p1 );
				if( kk == 0 )
					p1 = new Ycb;
			}
			p2 = p2->next;
		}
	}
	delete p1;			// delete the newnode create at last circle

	for( i = 2 ; i <= 2 * node_sum ; i++ )			//计算Li1 。L阵主对角为1，从第二项开始
	{
		p2 = l[i];
		if (p2 != NULL)
		{
			if (p2->lnxt == 1)
				p2->zh = p2->zh / ud[1];
		}
	}

	/***************************************** LU分解 **********************************/

	/*************************** L阵对角全为 1 **********************************/
	/*************************** 行号和 ycb 阵相同 ***************************/
	for( r = 2 ; r <= 2 * node_sum - 1 ; r++ ) 		 //从第二行开始, 算到倒数第二行
	{

		/****************************** compute Uri *****************************/
		for( i = 1 ; i <= r ; i++ )					// 每次循环存 L的 t数组先清零
			t[i] = 0;								// 【看L阵可知每行 l 元素的个数为 行号-1 】

		p2 = l[r];
		while( p2 != NULL )							//取 Lrk, r是行号，按行排的
		{
			t[p2->lnxt] = p2->zh;					// t数组为L阵列号
			p2 = p2->next;							// t内容为L的值
		}

		/******************************* U阵是除去对角的矩阵 ******************************/
		/************* U阵是对角元素单独存储。U阵行号与ycb相同但第一行没有 *******************/
		for( i = r + 1 ; i <= 2 * node_sum ; i++ )			// r为对角，要分开运算，所以从r+1开始
		{
			x = 0;
			p2 = u[i];										//从第 r + 1 列 U 阵开始。r列左边项为 Urr。 另行处理
			j = 0;											//
			while(( p2 != NULL ) && ( j == 0 ))				//取 Uki , k=(1~r-1) ,按列排的
			{
				if( p2->row >= r)							//U阵行号不可能大于列号，否则退出循环
					j = 1;
				else
				{
					x = x + t[p2->row] * p2->zh;			// Lrk*Uki ， u[i]的行对应l的列
					p2 = p2->next;
				}
			}

			if( x != 0 )
			{
				p1 = new Ycb;
				p1->row = r;
				p1->lnxt = i;
				p1->zh = - x;							//累加的和，求 Uri 的前期准备
				p1->next = NULL;
				u[p1->lnxt] = insert2( u[p1->lnxt] , p1 );
				if( kk == 1 )
					delete p1;
			}
		}

		/******************************* compute ud *****************************/
		x = 0;
		p2 = u[r];								//第r列 算对角u
		j = 0;
		while(( p2 != NULL ) && ( j == 0 ))		//取 Uki , k=(1~r-1) ,按列排的
		{
			if( p2->row >= r )					//U阵行号不可能大于列号，否则退出循环
				j = 1;
			else
			{
				x = x + t[p2->row] * p2->zh;	// Lrk * Uki
				p2 = p2->next;
			}
		}
		ud[r] = ud[r] - x;

		/******************************* compute L ****************************/

		for( i = 1 ; i <= r ; i++ )
			t[i] = 0;

		p2 = u[r];
		while( p2 != NULL )						//取Ukr,按列排的
		{
			t[p2->row] = p2->zh;
			p2 = p2->next;
		}
		for( i = r + 1 ; i <= 2 * node_sum ; i++ )
		{
			x = 0;
			p2 = l[i];
			j = 0;
			while(( p2 != NULL ) && ( j == 0 ))			//取Lik,k=(1~r-1)，按行排的
			{
				if( p2->lnxt >= r )
					j = 1;
				else
				{
					x = x + t[p2->lnxt] * p2->zh;		// Lik*Ukr , 最后统一除Urr
					p2 = p2->next;
				}
			}

			if( x != 0 )
			{
				p1 = new Ycb;
				p1->row = i;
				p1->lnxt = r;
				p1->zh = - x;
				p1->next = NULL;
				l[p1->row] = insert3( l[p1->row] , p1 );
				if(kk == 1)
					delete p1;
			}
		}
		for( i = r + 1 ; i <= 2 * node_sum ; i++ )			//   /Urr
		{
			p2 = l[i];
			j = 0;
			while(( p2 != NULL ) && ( j == 0 ))
			{
				if( p2->lnxt == r )
					p2->zh = p2->zh / ud[r];

				if( p2->lnxt >= r )
					j = 1;

				else
					p2 = p2->next;
			}
		}
	}

	/*************************** compute Unn对角元素 ****************************/
	for( i = 1 ; i <= 2 * node_sum ; i++ )
		t[i] = 0;

	p2 = l[ 2 * node_sum ];					//上面 i = 2*node_sum - 1 没用到的 l 最后一行
	while( p2 != NULL )						//取Lnk
	{
		t[p2->lnxt] = p2->zh;
			p2 = p2->next;
	}

	x = 0;
	p2 = u[ 2 * node_sum ];
	j = 0;
	while(( p2 != NULL ) && ( j == 0))		//取Uki,k=(1~r-1)
	{
		if( p2->row >= 2 * node_sum )
			j = 1;
		else
		{
			x = x + t[p2->row] * p2->zh;		//Lrk*Uki
			p2 = p2->next;
		}
	}
	ud[ 2 * node_sum ] = ud[ 2 * node_sum ] - x;	//ud[]原先存有 ycb 的最后一个对角元素

	/*************************** solvs Ly=b b=radv[] ***************************/
	for( i = 2 ; i <= 2 * node_sum ; i++ )
	{
		x = 0;
		j = 0;
		p2 = l[i];
		while( p2 != NULL )
		{
			x = x + radv[p2->lnxt] * p2->zh;	//Lik*yk,Lik是按行排列的
			p2 = p2->next;
		}
		radv[i] = radv[i] - x;
	}

	/****************************** solvs Ux=y y=radv[] ************************/
	radv[ 2 * node_sum ] = radv[ 2 * node_sum ] / ud[ 2 * node_sum ];
	for( i = 2 * node_sum - 1 ; i >= 1 ; i-- )
	{
		for( r = i + 1 ; r <= 2 * node_sum ; r++ )
			t[r] = 0;

		for( r = i + 1 ; r <= 2 * node_sum ; r++ )	//取Uir,这是按行排列的，要周转一下
		{
			p2 = u[r];
			j = 0;
			while(( p2 != NULL ) && ( j == 0 ))
			{
				if( p2->row == i )
					t[r] = p2->zh;
				if( p2->row >= i )
					j = 1;
				else
					p2 = p2->next;
			}
		}

		x = 0;
		for( r = i + 1 ; r <= 2 * node_sum ; r++ )
			x = x + t[r] * radv[r];
		radv[i] = (radv[i] - x) / ud[i];
	}

	for( i = 1 ; i <= 2 * node_sum ; i++ )			//释放空间
	{
		while( l[i] != NULL )
		{
			p2 = l[i]->next;
			delete l[i];
			l[i] = p2;				//l矩阵是按行排列的，其行号>列号
		}
	}
	for( i = 1 ; i <= 2 * node_sum ; i++ )		//释放空间
	{
		while( u[i] != NULL )
		{
			p2 = u[i]->next;
			delete u[i];
			u[i] = p2;				//u矩阵是按列排列的，其列号>行号
		}
	}
	return ;
}

/***************************************** 求解每节点的 △Pi △Gi *********************************************************/
void deltapq()
{

	struct Yz *p1 = NULL;
	struct pvd *p2 = NULL;
	for(int ii = 1 ; ii <= node_sum ; ii++)		//  △Pi  △Qi
	{
		p1 = y[ii];				//注意y阵中没有对角元素
		radv[2 * ii - 1] = 0;
		radv[2 * ii] =  0;
		while( p1 != NULL )						//原先radv[]中存的是0
		{
			radv[2 * ii - 1] += v[ii] * p1->y * v[p1->lnxt] * cos( rad[p1->row] - rad[p1->lnxt] - p1->rad );
			radv[2 * ii] += v[ii] * p1->y * v[p1->lnxt] * sin( rad[p1->row] - rad[p1->lnxt] - p1->rad );
			p1 = p1->next;
		}
		//再减去对角元素
		radv[2 * ii - 1] +=  v[ii] * ygbd[ii].g * v[ii] * cos( - ygbd[ii].b );
		radv[2 * ii] +=  v[ii] * ygbd[ii].g * v[ii] * sin( - ygbd[ii].b );

		radv[2 * ii - 1] = p[ii] - radv[2 * ii - 1];
		radv[2 * ii] = q[ii] - radv[2 * ii];
	}
	radv[2 * phd - 1] = 0;		// phd, △Pi  △Qi = 0
	radv[2 * phd] = 0;

	//******************************************** pvd
	p2 = headpvd;				//pvd, △Qi = 0
	while (p2 != NULL)
	{
		radv[2 * p2->node] = 0;
		p2 = p2->next;
	}

	for (int ii = 1; ii <= 2 * node_sum; ii++)
	{
		if (fabs(radv[ii]) > err)						//登记最大误差
			err = fabs(radv[ii]);
	}
}

/***************************************** 修正后的各个点的 相角θ 和 幅值v **********************************************/
void fixradv()
{
	for(int ii = 1 ; ii <= node_sum ; ii++)		// delta Pi || delta Qi
	{
		rad[ii] = rad[ii] + radv[2 * ii - 1];
		v[ii] = v[ii] + radv[2 * ii];
	}
}

/***************************************** 修正后的 平衡点的 Pi Gi 、 PV点的 Gi ************************************************/
void fixphdpvd()
{
	double x = 0, r = 0;

	struct Yz *p1 = NULL;								//平衡点处理
	p1 = y[phd];

	while( p1 != NULL )
	{
		x = x + p1->y * v[p1->lnxt] * cos(rad[phd] - rad[p1->lnxt] - p1->rad);
		r = r + p1->y * v[p1->lnxt] * sin(rad[phd] - rad[p1->lnxt] - p1->rad);
		p1 = p1->next;
	}
	/******************是否要把对角元素单独考虑？**************************/
	x = x + ygbd[phd].g * cos(rad[phd] - rad[phd] - ygbd[phd].b);
	r = r + ygbd[phd].g * sin(rad[phd] - rad[phd] - ygbd[phd].b);
	p[phd] = x * v[phd];
	q[phd] = r * v[phd];

	//*****************************PV点处理，p2为 pvd 链表***************************

	struct pvd *p2 = headpvd;
	struct Yz *p3 = NULL;
	p3 = y[p2->node];
	int bz = 0;

	while (p2 != NULL)
	{
		r = 0;
		p3 = y[p2->node];		// p2为 pvd 链表, 取其节点号，得到导纳阵表头
		while (p3 != NULL)
		{
			r = r + p3->y * v[p3->lnxt] * sin(rad[p3->row] - rad[p3->lnxt] - p3->rad);
			p3 = p3->next;
		}

		// 无功越界
		/*
		if (r > gmax[p2->node])
		{
			if (p2 == headpvd)
				bz = 1;
			q[p2->node] = gmax[p2->node];
			pvyjbz[p2->node] = 1;
			deletepvd(p2->node);
		}
		else if (r < gmin[p2->node])
		{
			if (p2 == headpvd)
				bz = 1;
			q[p2->node] = gmin[p2->node];
			pvyjbz[p2->node] = 1;
			deletepvd(p2->node);
		}
		else
			q[p2->node] = r * v[p2->node];
		p2 = p2->next;
		if (bz == 1)
			p2 = headpvd;*/
	}
	
}

/*****************************************删除指定的PV节点*********************************************************/
void deletepvd(int node0)
{
	struct pvd *p0 = headpvd;
	struct pvd *p1 = p0;
	while (p0 != NULL && p0->node != node0 )
	{
		p1 = p0;
		p0 = p0->next;
	}
	p1->next = p0->next;
	if (p0->node = headpvd->node)
		headpvd = p0->next;
	delete p0;
}

/*****************************************节点导纳矩阵插入指针数据*********************************************************/
struct Yz *insert1(struct Yz *tp, struct Yz *z)	//节点导纳矩阵插入指针数据
{
	struct Yz *p0, *p111, *p112 = NULL;
	double r, r1, x, x1;
	kk = 0;
	p111 = tp;			//指向原链表
	p112 = p111;
	p0 = z;				//指向要插入的那个节点

	if (p0 == NULL)
		return tp;
	if (tp == NULL)
	{
		tp = p0;
		p0->next = NULL;
		return tp;
	}
	if (p0->lnxt < p111->lnxt)			//插入的节点在链表开头之前
	{
		tp = p0;
		p0->next = p111;
		return tp;
	}
	while ((p0->lnxt > p111->lnxt) && (p111->next != NULL))
	{
		p112 = p111;
		p111 = p111->next;
	}

	if (p0->lnxt == p111->lnxt)		//两点间有多条线路或变压器
	{
		r = p111->y * cos(p111->rad);
		x = p111->y * sin(p111->rad);
		r1 = p0->y * cos(p0->rad);
		x1 = p0->y * sin(p0->rad);
		r = r + r1;
		x = x + x1;
		x1 = sqrt(r * r + x * x);
		p111->y = x1;
		p111->rad = atan2(x, r);
		kk = 1;
		return tp;
	}
	if ((p111->next == NULL) && (p0->lnxt > p111->lnxt))
	{
		p111->next = p0;
		p0->next = NULL;
	}
	else
	{
		p112->next = p0;
		p0->next = p111;
	}
	return tp;
}

/*****************************************雅可比矩阵插入指针数据 (按行 row)***********************************************/

struct Ycb *insert3(struct Ycb *tp, struct Ycb *z)	//节点导纳矩阵插入指针数据
{
	struct Ycb *p0, *p111, *p112 = NULL;

	kk = 0;
	p111 = tp;				//指向原链表
	p112 = p111;
	p0 = z;					//指向要插入的那个节点

	if (p0 == NULL)
		return tp;
	if (tp == NULL)
	{
		tp = p0;
		p0->next = NULL;
		return tp;
	}
	if (p0->lnxt < p111->lnxt)			//插入的节点在链表开头之前
	{
		tp = p0;
		p0->next = p111;
		return tp;
	}
	while ((p0->lnxt > p111->lnxt) && (p111->next != NULL))			//遍历链表抵达要插位置，或到末尾停止
	{
		p112 = p111;
		p111 = p111->next;
	}
	if (p0->lnxt == p111->lnxt)										//两点重合
	{
		p111->zh += p0->zh;
		kk = 1;
		return tp;
	}
	if ((p111->next == NULL) && (p0->lnxt > p111->lnxt))			//要插的位置位于表尾之后
	{
		p111->next = p0;
		p0->next = NULL;
	}
	else															//插入位置位于遍历抵达的位置
	{
		p112->next = p0;
		p0->next = p111;
	}
	return tp;
}


/*****************************************雅可比矩阵插入指针数据 (按列 lnxt)***********************************************/
struct Ycb *insert2(struct Ycb *tp, struct Ycb *z)	//节点导纳矩阵插入指针数据
{
	struct Ycb *p0, *p111, *p112 = NULL;

	kk = 0;
	p111 = tp;				//指向原链表
	p112 = p111;
	p0 = z;					//指向要插入的那个节点

	if (p0 == NULL)
		return tp;
	if (tp == NULL)
	{
		tp = p0;
		p0->next = NULL;
		return tp;
	}
	if (p0->row < p111->row)			//插入的节点在链表开头之前
	{
		tp = p0;
		p0->next = p111;
		return tp;
	}
	while ((p0->row > p111->row) && (p111->next != NULL))			//遍历链表抵达要插位置，或到末尾停止
	{
		p112 = p111;
		p111 = p111->next;
	}
	if (p0->row == p111->row)							//两点重合
	{
		p111->zh += p0->zh;
		kk = 1;
		return tp;
	}
	if ((p111->next == NULL) && (p0->row > p111->row))			//要插的位置位于表尾之后
	{
		p111->next = p0;
		p0->next = NULL;
	}
	else															//插入位置位于遍历抵达的位置
	{
		p112->next = p0;
		p0->next = p111;
	}
	return tp;
}
