//IC 2020
//ALOCAÇÃO DE CHAVES DE MANOBRA PARA RESTAURAÇÃO

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <time.h>
#include <cmath>
#include <complex>
#include <vector>
#include <string>
#include <algorithm>
#include <iterator>
#include <tuple>

using namespace std;

#define num_run_alg 30 //numero em que o algoritimo deve ser compilado

// DADOS DO SISTEMA DE 136 BARRAS ---------------------------------------

//Dados para o arquivo txt
#define linha_dados 157 //numero de linhas da matriz de dados +1
#define coluna_dados 9 //numero de colunas da matriz de dados +1

//Caracteristicas -----------------------------------------------------

#define num_AL 9 //numero de alimentadores +1 por conta do cpp
#define estado_restaurativo_pu 0.93 //tensao minima no estado restaurativo

//Chave a cada quantos kVA?
#define parametroCH_kVA 800;

//Parametros Funcao objetivo
#define tempo_falha 4 //numero de horas que o sistema fica em estado restaurativo
#define tempo_isolacao 0.06 //tempo necessario para fazer as manobras em horas
#define taxa_falhas 0.065 //taxa de falhas por km no ano
#define custoKWh 0.12 // em real 0.53187 (cotação 2017 ANEEL - Elektro - Sudeste)

//Caracteristicas Fluxo de Potencia ------------------------------------

//valores base, ou, de referencia
float sref = 100 * pow(10, 6); //100MVA
float vref = 13800; //13.8kV
float zref = (vref * vref) / sref;
float iref = sref / vref;

//tensao inicial nos nós do sistema: vref * (valor para dar a tensao inicial)
float pu_inicial_alimentador = 1.05;
float tensao_inicial_nos = (vref) * pu_inicial_alimentador;

//critério de convergencia do fluxo de potencia
complex <float> criterio_conv = 1 * pow(10, -5);
float epsilon = abs(criterio_conv);

int max_interacao = 8;

// Caracteristicas dos Alimentadores e Subestacoes ---------------------

int alimentadores[num_AL] = { 0, 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007 };
float capSE[num_AL] = { 0, sref, sref, sref, sref, sref, sref, sref, sref };

//NOTAS:
/*
- chaves no estado aberto eh considerado 0;
- chaves fechadas é considerado 1;
- VND com duas vizinhanças de intensificação;
- RVNS com quatro vizinhanças para diversificação;


*/

//#######################################################################################

// CLASSES ---------------------------------------------------

class ParametrosSistema
{
public:

	int noi[linha_dados] = {}; //nó inicial
	int nof[linha_dados] = {}; //nó final

	float lt_r[linha_dados] = {}; //parte real da lt
	float lt_x[linha_dados] = {}; //parte imaginaria da lt

	float s_nofr[linha_dados] = {}; //potencia complexa do nof real
	float s_nofq[linha_dados] = {}; //potencia complexa do nof img

	int candidato_aloc[linha_dados] = {}; //candidato a alocação de chaves
	int estado_swt[linha_dados] = {}; //estado da chave

	float dist_no[linha_dados] = {}; //distancia entre nós

	float potencia_al[num_AL] = {}; //potencia de cada alimentador

	float total_ativa = 0;
	complex <float> total_complexa = complex <float>(float(0.0), float(0.0));

	complex <float> lt[linha_dados] = {}; //linha de transmissao entre nós
	complex <float> s_nof[linha_dados] = {}; //potencia complexa do nof

	complex <float> pu_lt[linha_dados] = {}; //linha de transmissao entre nós em pu
	complex <float> pu_s_nof[linha_dados] = {}; //potencia complexa do nof em pu

	void leitura_parametros();
	void somatorio_potencia();

}ps;

class FluxoPotencia
{
public:

	int contadorFXP = 0; //conta quantas vzs realizou o processo de fluxo de potencia

	int camadaAL[num_AL][linha_dados][linha_dados] = {};

	int conexao_predef[linha_dados][3] = {};

	complex <float> tensao_inicial = complex <float>(float(tensao_inicial_nos / vref), float(0)); // tensao complexa nos nós na 1 iteraçao do fluxo de potencia

	complex <float> corrente_pu[linha_dados] = {};
	complex <float> tensao_pu[linha_dados] = {};

	void valores_nominais_tensao();
	void fluxo_potencia();

	void camadas(int alimentador, int camadaalimentador[linha_dados][linha_dados]);
	void backward_sweep(int camadaAL[linha_dados][linha_dados]);
	void forward_sweep(int alimentador, int camada[linha_dados][linha_dados]);

}fxp;

class AlocacaoChaves
{
public:

	int numch_AL[num_AL] = {}; //numero de chaves por alimentador seguindo o criterio estipulado
	int numch_SIS = 0;
	int posicaochaves[num_AL][linha_dados] = {}; //vetor com as posicoes das chaves
	int adjacente_chaves[num_AL][linha_dados][linha_dados] = {}; //secoes de todas as chaves, inclusive do disjuntor do alimentador 
	int secoes_chaves[num_AL][linha_dados][linha_dados] = {}; //secoes das chaves
	int chi[num_AL][linha_dados] = {}; //barra inicial da chave
	int chf[num_AL][linha_dados] = {}; //barra final da chave
	int antchi[num_AL][linha_dados] = {}; //barra inicial da chave
	int antchf[num_AL][linha_dados] = {}; //barra final da chave

	void chaves_anteriores();
	void volta_chaves_anteriores();
	void criterio_numero_de_chaves();
	void secoes_alimentador();
	float calculo_funcao_objetivo();

	tuple <int, float> contagem_criterio(int camada[linha_dados][linha_dados]); //criterio para a contagem de quantas chaves alocar em cada alimentador do sistema teste
	void adjacentes(int posicao[linha_dados], int adj[linha_dados][linha_dados], int alimentador); //calcula os adjacentes das chaves e da secao do alimentador
	float FO(float potencia_secao, float comprimeto_secao, float ens_isolacao);
	float energia_suprida(int AL, int SEC_al);

}ac;

class GVNS
{
public:

	//variaveis de solucao
	float current_solution = 0.0;
	float incumbent_solution = 0.0;
	float vnd_current = 0.0;
	float vnd_incumbent = 0.0;

	int analisepos = 0.0;
	float fo_al[num_AL] = {};

	//contadores
	int q_rvns1 = 0;
	int q_rvns2 = 0;
	int q_rvns3 = 0;
	int q_rvns4 = 0;
	int q_vnd1 = 0;
	int q_vnd2 = 0;

	//funcoes
	void primeiraaloc();
	float v1_RVNS(); //sortear duas chaves e aplicar VND
	float v2_RVNS(); //escolher todas as chaves de um alimentador sorteado e aplicar VND
	float v3_RVNS(); //selecionar dois alimentadores e sortear um numero de chaves para serem realocadas em posicoes aleatórias dentro do alimentador e fazer o VND em todas as chaves de cada alimentador
	float v4_RVNS(); //selecionar todos os alimentadores e manter apenas uma chave de cada alimentador na posicao original, o restante deve ser realocado em posicoes aleatorias dentro do alimentador

	float VND(vector<int>chaves);
	float v1_VND(vector<int>chavesv1); //mover para adjacente
	float v2_VND(vector<int>chavesv2); //mover para adjacente do adjacente

	void sorteiochaves(int numch, int camada[linha_dados][linha_dados], int posicao_camada[linha_dados], int alimentador); //sorteio inicial das chaves

}gvns;

//------------------------------------------------------------

void ParametrosSistema::leitura_parametros()
{
	FILE* arquivo;

	if ((arquivo = fopen("dados136.txt", "r")) == NULL)
		return;

	for (int i = 1; i < linha_dados; i++)
	{
		fscanf(arquivo, "%d%d%f%f%f%f%d%d%f", &ps.noi[i], &ps.nof[i], &ps.lt_r[i], &ps.lt_x[i], &ps.s_nofr[i], &ps.s_nofq[i], &ps.candidato_aloc[i], &ps.estado_swt[i], &ps.dist_no[i]);
	}

	fclose(arquivo);

	//dados complexos

	for (int j = 1; j < linha_dados; j++)
	{
		//dos ramos 

		ps.lt[j].real(ps.lt_r[j]);
		ps.lt[j].imag(ps.lt_x[j]);

		//das potencias complexas

		ps.s_nof[j].real(ps.s_nofr[j] * 1000); //a multiplicação por mil é porque os dados estao em kW
		ps.s_nof[j].imag(ps.s_nofq[j] * 1000); //a multiplicação por mil é porque os dados estao em kVAr
	}

	//dados em pu

	for (int k = 1; k < linha_dados; k++)
	{
		//dos ramos 
		ps.pu_lt[k] = ps.lt[k] / zref;

		//das potencias complexas
		ps.pu_s_nof[k] = ps.s_nof[k] / sref;
	}
}

void ParametrosSistema::somatorio_potencia()
{
	//aqui se soma a potencia ativa total, complexa e ativa

	for (int i = 1; i < linha_dados; i++)
	{
		ps.total_ativa = ps.total_ativa + ps.s_nofr[i];
		ps.total_complexa = ps.total_complexa + ps.s_nof[i];
	}
	;
}

void FluxoPotencia::camadas(int alimentador, int camadaalimentador[linha_dados][linha_dados])
{
	//zerar camada 
	bool add = false;

	for (int i = 1; i < linha_dados; i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			camadaalimentador[i][j] = 0;
		}
	}

	//define a camada do alimentador 

	camadaalimentador[1][1] = alimentador;

	int x = 2;
	int y = 1;


	add = false;
	//montando matriz adjacente

	for (int i = 1; i < linha_dados; i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			if (j != 0)
			{
				for (int k = 1; k < linha_dados; k++)
				{
					add = false;

					if (ps.noi[k] == camadaalimentador[i][j])
					{
						if (ps.estado_swt[k] == 1)
						{
							for (int m = 1; m < linha_dados; m++)
							{
								for (int p = 1; p < linha_dados; p++)
								{
									if (ps.nof[k] == camadaalimentador[m][p])
									{
										add = true; //ja tem o elemento na matriz
									}
								}
							}

							if (add == false)
							{
								camadaalimentador[x][y] = ps.nof[k];
								y++;
							}

						}
					}
				}
				//

				for (int k = 1; k < linha_dados; k++)
				{
					add = false;

					if (ps.nof[k] == camadaalimentador[i][j])
					{
						if (ps.estado_swt[k] == 1)
						{
							for (int m = 1; m < linha_dados; m++)
							{
								for (int p = 1; p < linha_dados; p++)
								{
									if (ps.noi[k] == camadaalimentador[m][p])
									{
										add = true; //ja tem o elemento na matriz
									}
								}
							}

							if (add == false)
							{
								camadaalimentador[x][y] = ps.noi[k];
								y++;
							}
						}
					}
				}
			}
		}
		x++;
		y = 1;
	}
}

void FluxoPotencia::backward_sweep(int camadaAL_bs[linha_dados][linha_dados])
{
	//atribuição das correntes

	for (int i = 1; i < linha_dados; i++)
	{
		for (int k = 1; k < linha_dados; k++)
		{
			for (int j = 1; j < linha_dados; j++)
			{
				if (ps.nof[i] == camadaAL_bs[k][j] && camadaAL_bs[k][j] != 0 && ps.estado_swt[i] == 1)
				{
					fxp.corrente_pu[i] = conj(ps.pu_s_nof[i] / fxp.tensao_pu[i]);
				}
			}
		}
	}

	// somatorio das correntes nos ramos

	for (int i = linha_dados; i > 0; i--)
	{
		for (int k = linha_dados; k > 0; k--)
		{
			if (camadaAL_bs[i][k] != 0)
			{
				for (int j = 1; j < linha_dados; j++)
				{
					if (camadaAL_bs[i][k] == ps.nof[j] && ps.estado_swt[j] == 1)
					{
						for (int o = 1; o < linha_dados; o++)
						{
							if (ps.nof[o] == ps.noi[j] && ps.estado_swt[o] == 1)
							{
								fxp.corrente_pu[o] += fxp.corrente_pu[j];
							}
						}
					}
				}
			}
		}
	}
}

void FluxoPotencia::forward_sweep(int alimentador, int camada[linha_dados][linha_dados])
{

	complex <float> unit = fxp.tensao_inicial; //complexo 

	//atribuindo a tensao para o restante dos nós

	for (int i = 1; i < linha_dados; i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			for (int k = 1; k < linha_dados; k++)
			{
				if (camada[i][j] == ps.nof[k] && ps.noi[k] != alimentador && camada[i][j] != 0 && ps.estado_swt[k] == 1)
				{
					for (int l = 1; l < linha_dados; l++)
					{
						if (ps.nof[l] == ps.noi[k] && ps.estado_swt[l] == 1)
						{
							fxp.tensao_pu[k] = fxp.tensao_pu[l] - (fxp.corrente_pu[k] * ps.pu_lt[k]);
						}
					}

				}

				else if (camada[i][j] == ps.nof[k] && ps.noi[k] == alimentador && camada[i][j] != 0 && ps.estado_swt[k] == 1)
				{
					fxp.tensao_pu[k] = unit - (fxp.corrente_pu[k] * ps.pu_lt[k]);
				}
			}
		}
	}
}

void FluxoPotencia::valores_nominais_tensao()
{
	//PARA ANALISAR O FLUXO

	cout << "\n\n";
	cout << "Fluxo de Potencia \n";
	cout << "No:" << "\t";
	cout << "V:" << "\t";
	cout << "ang:" << "\n";

	/*
	cout << "\t";
	cout << "I:" << "\t";
	cout << "ang:" << "\t";
	cout << "\n";
	*/

	for (int i = 1; i < linha_dados; i++)
	{

		cout << ps.nof[i] << "\t";
		cout << abs(fxp.tensao_pu[i] * vref) << "\t";
		cout << (arg(fxp.tensao_pu[i])) * 180 / 3.141592 << "\n";

		/*
		cout << abs(fxp.corrente_pu[i] * iref) << "\t";
		cout << arg(fxp.corrente_pu[i]) * 180 / 3.141592 << "\n";
		*/
	}
}

void FluxoPotencia::fluxo_potencia() //alterar conforme o numero de alimentadores, modificando quantas vezes cada funcao eh chamada
{
	bool criterio_satisfeito = false;

	complex <float> tensao_aux[linha_dados];
	complex <float> corrente_aux[linha_dados];

	float convergencia_tensao[linha_dados];
	float convergencia_corrente[linha_dados];

	int som = 0;
	int iteracao = 0;

	//zerando vetores
	for (int i = 0; i < linha_dados; i++)
	{
		tensao_aux[i] = 0.0;
		corrente_aux[i] = 0.0;
		fxp.tensao_pu[i] = 0.0;
		fxp.corrente_pu[i]  = 0.0;
	}
	//definir as camadas de fxp.camadas de cada alimentador

	for (int i = 1; i < num_AL; i++)
	{
		camadas(alimentadores[i], fxp.camadaAL[i]);
	}

	//tensoes iniciais nas barras

	for (int i = 1; i < linha_dados; i++)
	{
		if (ps.estado_swt[i] == 1)
		{
			tensao_pu[i] = fxp.tensao_inicial;
		}
		else
		{
			tensao_pu[i] = 0.0;
		}
	}

	//correntes iniciais nas barras

	for (int i = 1; i < linha_dados; i++)
	{
		corrente_pu[i] = 0.0;
	}

	///////////////////

	//FLUXO DE POTENCIA

	iteracao = 0;

	while (criterio_satisfeito == false)
	{
		som = 0;
		iteracao += 1;

		//copiando valores anteriores

		for (int i = 1; i < linha_dados; i++)
		{
			tensao_aux[i] = fxp.tensao_pu[i];
			corrente_aux[i] = fxp.corrente_pu[i];
		}

		//1 passo: BACKWARD

		for (int i = 1; i < num_AL; i++)
		{
			backward_sweep(fxp.camadaAL[i]);
		}


		//2 passo: FORWARD

		for (int i = 1; i < num_AL; i++)
		{
			forward_sweep(alimentadores[i], fxp.camadaAL[i]);
		}

		//comparação de critério satisfeito

		for (int k = 1; k < linha_dados; k++)
		{
			convergencia_tensao[k] = abs(fxp.tensao_pu[k]) - abs(tensao_aux[k]);
			convergencia_corrente[k] = abs(fxp.corrente_pu[k]) - abs(corrente_aux[k]);
		}

		for (int k = 1; k < linha_dados; k++)
		{
			if (abs(convergencia_tensao[k]) >= epsilon)
			{
				som += 1;
			}

			if (abs(convergencia_corrente[k]) >= epsilon)
			{
				som += 1;
			}
		}

		if (som == 0 || iteracao == max_interacao)
		{
			criterio_satisfeito = true;
		}
	}


	fxp.contadorFXP++; //conta o numero de vzs do processo do fluxo de potencia

}

tuple <int, float> AlocacaoChaves::contagem_criterio(int camada[linha_dados][linha_dados])
{
	int num_crit = 0;
	float num = 0.0;
	float potencia = 0.0;
	complex <float> potencia_S = 0.0;
	float modS = 0.0;

	//somando potencia
	for (int i = 1; i < linha_dados; i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			//segundo passo
			for (int k = 1; k < linha_dados; k++)
			{
				if (camada[i][j] == ps.nof[k])
				{
					potencia_S += ps.s_nof[k];
					potencia += ps.s_nofr[k];
				}
			}
		}
	}

	modS = abs(potencia_S);
	modS = modS / 1000;

	num = modS / parametroCH_kVA;

	num_crit = round(num);

	if (num_crit < 2) { num_crit = 2; }

	//encontando o numero estipulado
	return make_tuple(num_crit, potencia);
}

void AlocacaoChaves::criterio_numero_de_chaves() //alterar conforme o numero de alimentadores, modificando quantas vezes cada funcao eh chamada
{
	for (int i = 1; i < num_AL; i++)
	{
		tie(ac.numch_AL[i], ps.potencia_al[i]) = contagem_criterio(fxp.camadaAL[i]);
	}
}

void AlocacaoChaves::adjacentes(int posicao[linha_dados], int adj[linha_dados][linha_dados], int alimentador)
{
	int contline = 0; //contador da linha
	int contadorch = 0;
	int aux = 0; //auxiliar nas barras

	contline = 0;

	for (int i = 1; i < linha_dados; i++)
	{
		contadorch = 1;

		for (int j = 1; j < linha_dados; j++)
		{
			if (posicao[i] == j)
			{
				contline++;
				adj[contline][contadorch] = ps.nof[j];

				aux = 1;

				while (aux != linha_dados)
				{
					//localizando adjacentes
					for (int k = 1; k < linha_dados; k++)
					{
						if (ps.noi[k] == adj[contline][aux] && ps.candidato_aloc[k] == 1)
						{
							contadorch++;
							adj[contline][contadorch] = ps.nof[k];
						}
					}

					aux++;
				}
			}
		}

	}


	//agora pega a secao do alimentador

	contadorch = 1; //coloca o contador na primeira posicao

	for (int j = 1; j < linha_dados; j++)
	{
		if (alimentador == ps.noi[j])
		{
			contline++;
			adj[contline][contadorch] = ps.nof[j];

			aux = 1;

			while (aux != linha_dados)
			{
				//localizando adjacentes
				for (int k = 1; k < linha_dados; k++)
				{
					if (ps.noi[k] == adj[contline][aux] && ps.candidato_aloc[k] == 1)
					{
						contadorch++;
						adj[contline][contadorch] = ps.nof[k];
					}
				}

				aux++;
			}
		}
	}
}

void AlocacaoChaves::secoes_alimentador()
{
	int cont_ch[num_AL][linha_dados]; //contador de barras das chaves
	int cont = 0; //contador
	int pos = 0;

	//zerar matrizes dos adjacentes
	for (int i = 1; i < num_AL; i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			for (int k = 1; k < linha_dados; k++)
			{
				ac.secoes_chaves[i][j][k] = 0;
				ac.adjacente_chaves[i][j][k] = 0;
			}
		}
	}

	//encontrando todos os adjacentes
	for (int i = 1; i < num_AL; i++)
	{
		adjacentes(ac.posicaochaves[i], ac.adjacente_chaves[i], alimentadores[i]);
	}

	//zera contadores 
	for (int i = 1; i < num_AL; i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			cont_ch[i][j] = 0;
		}
	}

	//

	cont = 0;

	//determina qual das chaves tem mais barras adjacentes
	for (int i = 1; i < num_AL; i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			for (int k = 1; k < linha_dados; k++)
			{
				if (ac.adjacente_chaves[i][j][k] != 0)
				{
					cont++;
				}
			}
			cont_ch[i][j] = cont;
			cont = 0;
		}
	}

	//coloca em ordem a matriz antes de zerar


	for (int i = 1; i < num_AL; i++)
	{
		int xx = 0;
		xx = 0;
	agn:

		//pegando a maior
		cont = 0;
		pos = 1;

		for (int j = 1; j < linha_dados; j++)
		{
			if (cont < cont_ch[i][j])
			{
				cont = cont_ch[i][j];
				pos = j;
			}
		}

		cont_ch[i][pos] = 0;
		xx++;

		if (cont != 0)
		{
			for (int j = 1; j < linha_dados; j++)
			{
				ac.secoes_chaves[i][xx][j] = ac.adjacente_chaves[i][pos][j];
			}
			goto agn;
		}

		//aproveitar para colocar as barras adjacentes da chave em ordem

		//zerar para nao ter problemas futuros
		for (int j = 1; j < linha_dados; j++)
		{
			for (int k = 1; k < linha_dados; k++)
			{
				ac.adjacente_chaves[i][j][k] = 0;
			}

		}

		//ordenando...
		for (int j = 1; j < linha_dados; j++)
		{
			for (int k = 1; k < linha_dados; k++)
			{
				ac.adjacente_chaves[i][j][k] = ac.secoes_chaves[i][j][k];
			}

		}

		//agora sim, zerar os repetidos da array ac.adjacente_chaves[][][]

		for (int j = 1; j < linha_dados; j++)
		{
			for (int k = 1; k < linha_dados; k++)
			{
				//
				if (ac.secoes_chaves[i][j][k] == 0)
				{
					continue;
				}

				for (int y = 1; y < linha_dados; y++)
				{
					for (int p = 1; p < linha_dados; p++)
					{
						if (ac.secoes_chaves[i][y][p] == 0)
						{
							continue;
						}
						if (ac.secoes_chaves[i][j][k] == ac.secoes_chaves[i][y][p] && y != j)
						{
							ac.secoes_chaves[i][j][k] = 0;
						}
					}
				}
			}
		}
	}
}

float AlocacaoChaves::energia_suprida(int AL, int SEC_al)
{
	bool analise;
	float energia_sup = 0.0;
	complex <float> capacidadeSE;

	

	// ANALISANDO FLUXO DE POTENCIA:
	
	//faz o fluxo de potencia para a situação do remanejamento
	fxp.fluxo_potencia();

	analise = true; //verdadeiro se atende as condiçoes, falso caso contrario - assume-se inicialmente

	// 1. Níveis de tensao devem estar acima de 0.93pu
	for (int y = 1; y < linha_dados; y++)
	{
		if (abs(fxp.tensao_pu[y]) < estado_restaurativo_pu && abs(fxp.tensao_pu[y]) != 0)
		{
			analise = false;
			break;
		}
	}

	// 2. A potencia alimentada por cada alimentador deve estar dentro da capacidade do alimentador
	for (int i = 1; i < num_AL; i++)
	{
		capacidadeSE.real(0.0);
		capacidadeSE.imag(0.0);

		for (int j = 1; j < linha_dados; j++)
		{
			for (int k = 1; k < linha_dados; k++)
			{
				for (int t = 1; t < linha_dados; t++)
				{
					if (fxp.camadaAL[i][j][k] == ps.nof[t] && fxp.camadaAL[i][j][k] != 0 && ps.estado_swt[t] == 1)
					{
						capacidadeSE += ps.s_nof[t];
					}
				}
			}
		}

		if (abs(capacidadeSE) > capSE[i])
		{
			analise = false;
			break;
		}
	}

	if (analise == true)
	{
		//a energia suprida é a da secao
		energia_sup = 0.0;

		for (int i = 1; i < linha_dados; i++)
		{
			for (int j = 1; j < linha_dados; j++)
			{
				if (ps.nof[i] == ac.secoes_chaves[AL][SEC_al][j])
				{
					energia_sup = energia_sup + ps.s_nofr[i];
				}
			}
		}

		//caso: energia_sup = 0.0 --> atribuir -1.0 a essa variavel afim de zera-la depois
		if (energia_sup == 0)
		{
			energia_sup = -1.0;
		}
		
	}
	else
	{
		energia_sup = 0.0;
	}

	//retornar ao estado inicial
	ps.leitura_parametros();
	fxp.fluxo_potencia();

	//retornar valor a variavel
	return(energia_sup);
}

float AlocacaoChaves::FO(float potencia_secao, float comprimento, float ens_isolacao)
{
	float resultado = 0.0;

	if (ens_isolacao == potencia_secao)
	{
		//neste caso, nao adianta fazer manobras, a secao ficará desligada durante todo o processo de falta

		resultado = taxa_falhas * comprimento * (ens_isolacao * custoKWh * tempo_falha);
	}
	else
	{
		//neste caso, adianta o chaveamento, entao, terá o remanejamento de cargas

		resultado = (taxa_falhas * comprimento * (potencia_secao * custoKWh * tempo_falha)) + (ens_isolacao * tempo_isolacao);
	}


	return(resultado);
}

float AlocacaoChaves::calculo_funcao_objetivo()
{
	float comprimento_secao = 0.0;
	float potencia_W = 0.0;
	float valorFO = 0.0;
	float chamadaFO = 0.0;
	float potencia_isolacao = 0.0;
	float ENSotima = 0.0;
	float ens = 0.0;
	
	float valorFOtot = 0.0;

	int contcondicao = 0;

	vector <int> secao;
	vector <int> posicao;
	vector <vector<int>> analise_remanejamento;
	vector <vector<int>> remanej_cargas;

	posicao.clear();
	secao.clear();
	analise_remanejamento.clear();
	remanej_cargas.clear();

	//deve-se analisar todas as secoes para os valores da funcao obj
	ac.secoes_alimentador();

	//alimentador i
	for (int i = 1; i < num_AL; i++)
	{
		if (gvns.analisepos != i) { continue; }

		valorFO = 0.0;

		//secao j
		for (int j = 1; j < linha_dados; j++)
		{
			//ver se vale a pena fazer o laço, se o vetor estiver zerado é só custo computacional a toa
			contcondicao = 0;

			for (int k = 0; k < linha_dados; k++)
			{
				if (ac.secoes_chaves[i][j][k] != 0)
				{
					contcondicao++;
				}
			}

			if (contcondicao == 0) { continue; }

			//////////////////////

			comprimento_secao = 0.0;
			potencia_W = 0.0;
			potencia_isolacao = 0.0;
			ENSotima = 0.0;

			// k = barras da seção j

			//analise comprimento e potencia nao suprida
			for (int k = 1; k < linha_dados; k++)
			{
				//comprimento
				for (int y = 1; y < linha_dados; y++)
				{
					if (ac.secoes_chaves[i][j][k] == ps.nof[y])
					{
						comprimento_secao = comprimento_secao + ps.dist_no[y];
						ENSotima = ENSotima + ps.s_nofr[y];
					}
				}
			}

			// 1) primeiro deve-se pegar toda a area do alimentador e deliga-la
			for (int k = 1; k < linha_dados; k++)
			{
				for (int y = 1; y < linha_dados; y++)
				{
					if (ac.adjacente_chaves[i][1][k] == ps.nof[y])
					{
						potencia_isolacao = potencia_isolacao + ps.s_nofr[y];
					}
				}
			}


			// 2) agora deve-se fazer o devido chaveamento
			// 2a) isolando secao j
			for (int k = 1; k < linha_dados; k++)
			{
				for (int y = 1; y < linha_dados; y++)
				{
					if (ac.secoes_chaves[i][j][k] == ps.nof[y] && ac.secoes_chaves[i][j][k] != 0)
					{
						ps.estado_swt[y] = 0;
					}
				}
			}

			/*
			// 2b) definindo chaves para remanejamento
			posicao.clear();
			remanej_cargas.clear();
			bool adiciona_vetpos = true;

			for (int k = 1; k < linha_dados; k++)
			{
				if (k <= j) { continue; }

				for (int t = 0; t < linha_dados; t++)
				{
					if (ac.adjacente_chaves[i][k][t] == 0) { continue; }

					for (int n = 1; n < linha_dados; n++)
					{
						if (ps.nof[n] == ac.adjacente_chaves[i][k][t] && ps.candidato_aloc[n] == 0 && ps.estado_swt[n] == 0 && ps.noi[n] != alimentadores[i] && ps.nof[n] != alimentadores[i])
						{
							posicao.push_back(n);
						}

						else if (ps.noi[n] == ac.adjacente_chaves[i][k][t] && ps.candidato_aloc[n] == 0 && ps.estado_swt[n] == 0 && ps.noi[n] != alimentadores[i] && ps.nof[n] != alimentadores[i])
						{
							posicao.push_back(n);
						}
						
					}
				}

				if (posicao.size() != 0)
				{
					if (remanej_cargas.size() == 0)
					{
						//se estiver vazio, adiciona o que tiver encontrado
						remanej_cargas.push_back(posicao);
						posicao.clear();
					}
					else
					{
						//caso nao tiver tamanho nulo, analisar se ja tem esta opcao - se sim, nao adicciona - se não, adiciona em outra coluna
						adiciona_vetpos = true;
						for (int m = 0; m < remanej_cargas.size(); m++)
						{
							for (int v = 0; v < remanej_cargas[m].size(); v++)
							{
								for (int z = 0; z < posicao.size(); z++)
								{
									if(remanej_cargas[m][v] == posicao[z])
									{
										adiciona_vetpos = false; //ja tem esta posicao, nao add
									}
								}
							}
						}

						if (adiciona_vetpos == true)
						{
							//adiciona
							remanej_cargas.push_back(posicao);
							posicao.clear();
						}
						else
						{
							posicao.clear();
						}

					}

				}
			}
			
			*/
			
			// 2b) cenario da falta

			//zerar falha na camada
			for (int k = 1; k < linha_dados; k++)
			{
				for (int y = 1; y < linha_dados; y++)
				{
					for (int t = 1; t < linha_dados; t++)
					{
						if (fxp.camadaAL[i][k][y] == ac.secoes_chaves[i][j][t])
						{
							fxp.camadaAL[i][k][y] = 0;
						}
					}
				}
			}

			bool inicio;
			int cont_nao0;

			inicio = false;

			for (int k = 1; k < linha_dados; k++)
			{
				cont_nao0 = 0;

				for (int y = 1; y < linha_dados; y++)
				{
					if (fxp.camadaAL[i][k][y] != 0)
					{
						secao.push_back(fxp.camadaAL[i][k][y]);
						cont_nao0++;
					}
				}

				if (cont_nao0 == 0)
				{
					secao.clear();
					inicio = true;
				}
				else if (inicio == true && secao.size() != 0)
				{
					for (int m = 1; m < linha_dados; m++)
					{
						for (int n = 1; n < linha_dados; n++)
						{
							for (int t = 0; t < secao.size(); t++)
							{
								if (secao[t] == ac.secoes_chaves[i][m][n])
								{
									posicao.push_back(m);
								}
							}
						}
					}

					secao.clear();
				}
				else
				{
					secao.clear();
				}
			}

			//eliminar elementos iguais no vetor
			for (int k = 0; k < posicao.size(); k++)
			{
				for (int t = 0; t < posicao.size(); t++)
				{
					if (t != k && posicao[k] == posicao[t])
					{
						posicao[t] = 0;
					}
				}
			}

			//analisando o remanejamento:
			secao.clear();

			for (int k = 0; k < posicao.size(); k++)
			{
				if (posicao[k] == 0) { continue; }

				for (int t = 1; t < linha_dados; t++)
				{
					if (ac.adjacente_chaves[i][posicao[k]][t] != 0)
					{
						secao.push_back(ac.adjacente_chaves[i][posicao[k]][t]);
					}
				}

				//verificar se esta contido na camada
				bool contbar;

				contbar = false;

				for (int g = 0; g < analise_remanejamento.size(); g++)
				{
					for (int h = 0; h < analise_remanejamento[g].size(); h++)
					{
						for (int t = 0; t < secao.size(); t++)
						{
							if (analise_remanejamento[g][h] == secao[t]) { contbar = true; }
						}
					}
				}

				if (contbar == false)
				{
					analise_remanejamento.push_back(secao);
					secao.clear();
				}
				else
				{
					secao.clear();
				}


			}

			//pegando posições
			posicao.clear();

			for (int k = 0; k < analise_remanejamento.size(); k++)
			{
				for (int t = 0; t < analise_remanejamento[k].size(); t++)
				{
					for (int y = 1; y < linha_dados; y++)
					{
						if (ps.noi[y] == analise_remanejamento[k][t] && ps.candidato_aloc[y] == 0)
						{
							posicao.push_back(y);
						}
						else if (ps.nof[y] == analise_remanejamento[k][t] && ps.candidato_aloc[y] == 0)
						{
							posicao.push_back(y);
						}
					}
				}

				if (!posicao.empty())
				{
					remanej_cargas.push_back(posicao);
					posicao.clear();
				}
			}

			//3) Calculo da ENS pelo sistema caso ocorra falha na seção j do alimentador i

			if (remanej_cargas.size() != 0)
			{
				ens = ps.potencia_al[i];

				//somente uma opcao para remanejamento de cada vez
				for (auto& linha : remanej_cargas)
				{
					for (auto& coluna : linha)
					{
						ps.estado_swt[coluna] = 1;

						potencia_W = ac.energia_suprida(i, j);

						ps.estado_swt[coluna] = 0;

						if (potencia_W != 0)
						{
							break;
						}
					}

					if (potencia_W < 0.0)
					{
						potencia_W = 0.0; //do desvio da funcao energia_suprida
					}

					ens = ens - potencia_W;
					potencia_W = 0.0;
				}
			}
			else
			{
				//nao tem como fazer manobra, a ENS será os adjacentes da chave
				ens = 0.0;

				for (int y = 1; y < linha_dados; y++)
				{
					for (int h = 1; h < linha_dados; h++)
					{
						if (ps.nof[h] == ac.adjacente_chaves[i][j][y])
						{
							ens = ens + ps.s_nofr[h];
						}
					}
				}
			}

			analise_remanejamento.clear();
			remanej_cargas.clear();
			secao.clear();
			posicao.clear();


			// 4) chamando a funcao objetivo
			chamadaFO = 0.0;

			chamadaFO = ac.FO(ens, comprimento_secao, potencia_isolacao);

			valorFO += chamadaFO;
		}

		gvns.fo_al[i] = valorFO;
	}

	//soma
	valorFOtot = 0.0;
	for (int r = 1; r < num_AL; r++)
	{
		valorFOtot = valorFOtot + gvns.fo_al[r];
	}

	//imprime valor caso quiser
	//cout << "FO: " << valorFOtot << endl;

	return(valorFOtot);
}

void AlocacaoChaves::chaves_anteriores()
{
	for (int i = 1; i < num_AL; i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			ac.antchi[i][j] = ac.chi[i][j];
			ac.antchf[i][j] = ac.chf[i][j];
		}
	}
}

void AlocacaoChaves::volta_chaves_anteriores()
{
	for (int i = 1; i < num_AL; i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			ac.chi[i][j] = ac.antchi[i][j];
			ac.chf[i][j] = ac.antchf[i][j];
		}
	}
}

void GVNS::sorteiochaves(int numch, int camada[linha_dados][linha_dados], int posicao_camada[linha_dados], int alimentador)
{
	int barras_camada[linha_dados];
	int aux = 0;
	int sort = 0;
	bool atribuir = false; //variavel auxiliar
	bool igual = false;

	//zera barra_camada
	for (int i = 1; i < linha_dados; i++)
	{
		barras_camada[i] = 0;
	}

	//atribui as barras do alimentador
	aux = 1;

	for (int i = 1; i < linha_dados; i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			if (camada[i][j] != 0)
			{
				barras_camada[aux] = camada[i][j];
				aux++;
			}
		}
	}

	//sorteia
sorteio:
	for (int i = 1; i <= numch; i++)
	{
	rand_dnv:

		atribuir = false;

		sort = rand() % aux + 1;

		sort = barras_camada[sort];

		if (sort != 0)
		{
			//localizar posicao
			for (int j = 1; j < linha_dados; j++)
			{
				if (sort == ps.nof[j] && ps.candidato_aloc[j] == 1 && ps.noi[j] != alimentador)
				{
					posicao_camada[i] = j;
					atribuir = true;
				}
			}

			if (atribuir == false)
			{
				goto rand_dnv;
			}
		}
		else
		{
			goto rand_dnv;
		}
	}

	//analisa se as tres sao diferentes
	igual = false;
	for (int i = 1; i < linha_dados; i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			if (posicao_camada[i] == posicao_camada[j] && posicao_camada[i] != 0 && i != j)
			{
				igual = true;
			}
		}
	}

	if (igual == true)
	{
		goto sorteio;
	}
}

void GVNS::primeiraaloc() //alterar conforme o numero de alimentadores, modificando quantas vezes cada funcao eh chamada
{

	for (int i = 1; i < num_AL; i++)
	{
		sorteiochaves(ac.numch_AL[i], fxp.camadaAL[i], ac.posicaochaves[i], alimentadores[i]);
	}

}

float GVNS::v1_VND(vector <int> chavesv1)
{
	//troca as chaves para o vizinho, ou seja, para o adjacente
	//chavesv1 é o numero das chaves

	int identf = 0;
	vector <int> pos_vizinhas;

	bool repetepos = false;
	int pos_vetorcaso_nao_melhor = 0;

	float soluc = 0.0;
	float solucao = 0.0;

	vector<int>auxfunc; //vetor auxiliar das funcoes
	vector <vector<int>> identch; //sempre tera duas posicoes:: 1)chi - 2)chf

	solucao = gvns.vnd_incumbent;

	// 1) localizar chaves sorteadas

	for (int i = 0; i < chavesv1.size(); i++)
	{
		identf = 0;

		for (int j = 1; j < num_AL; j++)
		{
			for (int k = 1; k < linha_dados; k++)
			{
				if (ac.chi[j][k] != 0 && ac.chf[j][k] != 0)
				{
					identf++;

					if (identf == chavesv1[i])
					{
						auxfunc.push_back(ac.chi[j][k]);
						auxfunc.push_back(ac.chf[j][k]);

						identch.push_back(auxfunc);
						auxfunc.clear();
					}
				}
			}
		}
	}

	//2) identificar chaves vizinhas e colocalas em um vetor e executar demais passos

	for (int i = 0; i < identch.size(); i++)
	{
		for (int k = 1; k < linha_dados; k++)
		{
			if (identch[i][0] == ps.noi[k])
			{
				if (identch[i][1] == ps.nof[k]) { continue; }
				else
				{
					repetepos = false;
					for (int j = 0; j < pos_vizinhas.size(); j++)
					{
						if (pos_vizinhas[j] == k) { repetepos = true; }
					}

					if (repetepos == false)
					{
						if (ps.candidato_aloc[k] == 1)
						{
							pos_vizinhas.push_back(k);
						}
					}
				}
			}

			if (identch[i][0] == ps.nof[k])
			{
				if (identch[i][1] == ps.noi[k]) { continue; }
				else
				{
					repetepos = false;
					for (int j = 0; j < pos_vizinhas.size(); j++)
					{
						if (pos_vizinhas[j] == k) { repetepos = true; }
					}

					if (repetepos == false)
					{
						if (ps.candidato_aloc[k] == 1)
						{
							pos_vizinhas.push_back(k);
						}
					}

				}
			}

			if (identch[i][1] == ps.noi[k])
			{
				if (identch[i][0] == ps.nof[k]) { continue; }
				else
				{
					repetepos = false;
					for (int j = 0; j < pos_vizinhas.size(); j++)
					{
						if (pos_vizinhas[j] == k) { repetepos = true; }
					}

					if (repetepos == false)
					{
						if (ps.candidato_aloc[k] == 1)
						{
							pos_vizinhas.push_back(k);
						}
					}
				}
			}

			if (identch[i][1] == ps.nof[k])
			{
				if (identch[i][0] == ps.noi[k]) { continue; }
				else
				{
					repetepos = false;
					for (int j = 0; j < pos_vizinhas.size(); j++)
					{
						if (pos_vizinhas[j] == k) { repetepos = true; }
					}

					if (repetepos == false)
					{
						if (ps.candidato_aloc[k] == 1)
						{
							pos_vizinhas.push_back(k);
						}
					}
				}
			}
		}

		//3) Analisar qual o melhor vizinho

		int pos1ch = 0;
		int pos2ch = 0;

		//identificando chave para mudar
		for (int t = 1; t < num_AL; t++)
		{
			for (int a = 1; a < linha_dados; a++)
			{
				if (ac.chi[t][a] == identch[i][0] && ac.chf[t][a] == identch[i][1])
				{
					pos1ch = t;
					pos2ch = a;
				}
			}
		}

		//analisar a melhor posicao - IMPORTANTE -
		soluc = 0.0;
		pos_vetorcaso_nao_melhor = ac.posicaochaves[pos1ch][pos2ch];

		for (int t = 0; t < pos_vizinhas.size(); t++)
		{
			ac.chaves_anteriores(); //salva chaves anteriores

			//atribuir nova chave
			ac.chi[pos1ch][pos2ch] = ps.noi[pos_vizinhas[t]];
			ac.chf[pos1ch][pos2ch] = ps.nof[pos_vizinhas[t]];

			//atribuir nova posicao da chave
			ac.posicaochaves[pos1ch][pos2ch] = pos_vizinhas[t];

			//refazer novas camadas
			ac.secoes_alimentador();

			gvns.analisepos = pos1ch;
			soluc = ac.calculo_funcao_objetivo();

			//analisar caso
			if (soluc > solucao)
			{
				//volta oq estava
				ac.volta_chaves_anteriores();
				ac.posicaochaves[pos1ch][pos2ch] = pos_vetorcaso_nao_melhor;
				ac.secoes_alimentador();
			}
			else
			{
				//atualiza
				solucao = soluc;
			}
		}

		pos_vizinhas.clear();

	}
	
	chavesv1.clear();
	
	return solucao;
}

float GVNS::v2_VND(vector <int> chavesv2)
{
	//troca as chaves para o vizinho do vizinho, ou seja, adjacente do adjacente
	//chavesv2 é o numero das chaves

	int identf = 0;
	vector <int> pos_vizinhas;

	bool repetepos = false;
	int pos_vetorcaso_nao_melhor = 0;

	vector <float> result_parcial; //resultado parcial da funcao objetivo
	float soluc = 0.0;
	float solution = 0.0;

	vector<int>auxfunc; //vetor auxiliar das funcoes
	vector <vector<int>> identch; //sempre tera duas posicoes:: 1)chi - 2)chf
	vector <vector<int>> identch_aux; //sempre tera duas posicoes:: 1)chi - 2)chf - eh o vetor auxiliar usado para encontrar as chaves

	solution = gvns.vnd_incumbent;
	// 1) localizar chaves sorteadas

	for (int i = 0; i < chavesv2.size(); i++)
	{
		identf = 0;

		for (int j = 1; j < num_AL; j++)
		{
			for (int k = 1; k < linha_dados; k++)
			{
				if (ac.chi[j][k] != 0 && ac.chf[j][k] != 0)
				{
					identf++;

					if (identf == chavesv2[i])
					{
						auxfunc.push_back(ac.chi[j][k]);
						auxfunc.push_back(ac.chf[j][k]);

						identch.push_back(auxfunc);
						auxfunc.clear();
					}
				}
			}
		}
	}

	//2) identificar chaves vizinhas e colocalas em um vetor e executar demais passos

	for (int i = 0; i < identch.size(); i++)
	{
		//faz o adjacente

		for (int k = 1; k < linha_dados; k++)
		{
			if (identch[i][0] == ps.noi[k])
			{
				if (identch[i][1] == ps.nof[k]) { continue; }
				else
				{
					repetepos = false;
					for (int j = 0; j < pos_vizinhas.size(); j++)
					{
						if (pos_vizinhas[j] == k) { repetepos = true; }
					}

					if (repetepos == false)
					{
						if (ps.candidato_aloc[k] == 1)
						{
							pos_vizinhas.push_back(k);
						}
					}
				}
			}

			if (identch[i][0] == ps.nof[k])
			{
				if (identch[i][1] == ps.noi[k]) { continue; }
				else
				{
					repetepos = false;
					for (int j = 0; j < pos_vizinhas.size(); j++)
					{
						if (pos_vizinhas[j] == k) { repetepos = true; }
					}

					if (repetepos == false)
					{
						if (ps.candidato_aloc[k] == 1)
						{
							pos_vizinhas.push_back(k);
						}
					}

				}
			}

			if (identch[i][1] == ps.noi[k])
			{
				if (identch[i][0] == ps.nof[k]) { continue; }
				else
				{
					repetepos = false;
					for (int j = 0; j < pos_vizinhas.size(); j++)
					{
						if (pos_vizinhas[j] == k) { repetepos = true; }
					}

					if (repetepos == false)
					{
						if (ps.candidato_aloc[k] == 1)
						{
							pos_vizinhas.push_back(k);
						}
					}
				}
			}

			if (identch[i][1] == ps.nof[k])
			{
				if (identch[i][0] == ps.noi[k]) { continue; }
				else
				{
					repetepos = false;
					for (int j = 0; j < pos_vizinhas.size(); j++)
					{
						if (pos_vizinhas[j] == k) { repetepos = true; }
					}

					if (repetepos == false)
					{
						if (ps.candidato_aloc[k] == 1)
						{
							pos_vizinhas.push_back(k);
						}
					}
				}
			}
		}

		//atribuir as chaves e realizar novamente o processo
		auxfunc.clear();

		for (int t = 0; t < pos_vizinhas.size(); t++)
		{
			auxfunc.push_back(ps.noi[pos_vizinhas[t]]);
			auxfunc.push_back(ps.nof[pos_vizinhas[t]]);

			identch_aux.push_back(auxfunc);
			auxfunc.clear();
		}

		pos_vizinhas.clear(); //limpa as posicoes para pegar as novas

		// faz o adjacente dos adjacentes

		//copiar passo de identificação para pegar os proximos adjacentes
		bool partedachave = false;

		for (int t = 0; t < identch_aux.size(); t++)
		{
			for (int k = 1; k < linha_dados; k++)
			{
				if (identch_aux[t][0] == ps.noi[k])
				{
					if (identch_aux[t][1] == ps.nof[k]) { continue; }
					else
					{
						partedachave = false;
						for (int j = 0; j < identch[i].size(); j++)
						{
							if (identch[i][j] == ps.noi[k] || identch[i][j] == ps.nof[k])
							{
								partedachave = true;
							}
						}

						if (partedachave == false)
						{
							repetepos = false;
							for (int j = 0; j < pos_vizinhas.size(); j++)
							{
								if (pos_vizinhas[j] == k) { repetepos = true; }
							}

							if (repetepos == false)
							{
								if (ps.candidato_aloc[k] == 1)
								{
									pos_vizinhas.push_back(k);
								}
							}
						}
					}
				}

				if (identch_aux[t][0] == ps.nof[k])
				{
					if (identch_aux[t][1] == ps.noi[k]) { continue; }
					else
					{
						partedachave = false;
						for (int j = 0; j < identch[i].size(); j++)
						{
							if (identch[i][j] == ps.noi[k] || identch[i][j] == ps.nof[k])
							{
								partedachave = true;
							}
						}

						if (partedachave == false)
						{
							repetepos = false;
							for (int j = 0; j < pos_vizinhas.size(); j++)
							{
								if (pos_vizinhas[j] == k) { repetepos = true; }
							}

							if (repetepos == false)
							{
								if (ps.candidato_aloc[k] == 1)
								{
									pos_vizinhas.push_back(k);
								}
							}
						}

					}
				}

				if (identch_aux[t][1] == ps.noi[k])
				{
					if (identch_aux[t][0] == ps.nof[k]) { continue; }
					else
					{
						partedachave = false;
						for (int j = 0; j < identch[i].size(); j++)
						{
							if (identch[i][j] == ps.noi[k] || identch[i][j] == ps.nof[k])
							{
								partedachave = true;
							}
						}

						if (partedachave == false)
						{
							repetepos = false;
							for (int j = 0; j < pos_vizinhas.size(); j++)
							{
								if (pos_vizinhas[j] == k) { repetepos = true; }
							}

							if (repetepos == false)
							{
								if (ps.candidato_aloc[k] == 1)
								{
									pos_vizinhas.push_back(k);
								}
							}
						}
					}
				}

				if (identch_aux[t][1] == ps.nof[k])
				{
					if (identch_aux[t][0] == ps.noi[k]) { continue; }
					else
					{
						partedachave = false;
						for (int j = 0; j < identch[i].size(); j++)
						{
							if (identch[i][j] == ps.noi[k] || identch[i][j] == ps.nof[k])
							{
								partedachave = true;
							}
						}

						if (partedachave == false)
						{
							repetepos = false;
							for (int j = 0; j < pos_vizinhas.size(); j++)
							{
								if (pos_vizinhas[j] == k) { repetepos = true; }
							}

							if (repetepos == false)
							{
								if (ps.candidato_aloc[k] == 1)
								{
									pos_vizinhas.push_back(k);
								}
							}
						}
					}
				}
			}
		}

		//3) Analisar qual o melhor vizinho

		int pos1ch = 0;
		int pos2ch = 0;

		//identificando chave para mudar
		for (int t = 1; t < num_AL; t++)
		{
			for (int a = 1; a < linha_dados; a++)
			{
				if (ac.chi[t][a] == identch[i][0] && ac.chf[t][a] == identch[i][1])
				{
					pos1ch = t;
					pos2ch = a;
				}
			}
		}

		//analisar a melhor posicao - IMPORTANTE -
		soluc = 0.0;
		pos_vetorcaso_nao_melhor = ac.posicaochaves[pos1ch][pos2ch];

		for (int t = 0; t < pos_vizinhas.size(); t++)
		{
			ac.chaves_anteriores(); //salva chaves anteriores

			//atribuir nova chave
			ac.chi[pos1ch][pos2ch] = ps.noi[pos_vizinhas[t]];
			ac.chf[pos1ch][pos2ch] = ps.nof[pos_vizinhas[t]];

			//atribuir nova posicao da chave
			ac.posicaochaves[pos1ch][pos2ch] = pos_vizinhas[t];

			//refazer novas camadas
			ac.secoes_alimentador();

			gvns.analisepos = pos1ch;
			soluc = ac.calculo_funcao_objetivo();

			//analisar caso
			if (soluc > solution)
			{
				ac.volta_chaves_anteriores();
				ac.posicaochaves[pos1ch][pos2ch] = pos_vetorcaso_nao_melhor;
				ac.secoes_alimentador();
			}
			else
			{
				solution = soluc;
			}
		}

		pos_vizinhas.clear();
		auxfunc.clear();
		identch_aux.clear();
	}

	chavesv2.clear();

	return solution;

}

float GVNS::VND(vector <int> chaves)
{
	//O VND tem o objetivo de intensificar a busca do algoritmo
	
	gvns.vnd_incumbent = gvns.incumbent_solution;
	gvns.vnd_current = gvns.vnd_incumbent;

	vector<int>vetaux;

	//for (int i = 0; i < chaves.size(); i++)
	//{
	
		//vetaux.push_back(chaves[i]);

	inicioVND:

		gvns.vnd_current = gvns.v1_VND(chaves);

		if (gvns.vnd_current < gvns.vnd_incumbent)
		{
			gvns.vnd_incumbent = gvns.vnd_current;
			gvns.q_vnd1++;

			cout << "_ 1-VND: " << gvns.vnd_incumbent << endl;

			goto inicioVND;
		}

		gvns.vnd_current = gvns.v2_VND(chaves);

		if (gvns.vnd_current < gvns.vnd_incumbent)
		{
			gvns.vnd_incumbent = gvns.vnd_current;
			gvns.q_vnd2++;

			cout << "_ 2-VND: " << gvns.vnd_incumbent << endl;

			goto inicioVND;
		}

		//vetaux.clear();
	//}


	chaves.clear();

	//fim vnd
	return gvns.vnd_incumbent;
}

float GVNS::v1_RVNS()
{
	//Descrição: sortear duas chaves quaisquer alocadas no sistema e aplicar o VND
	//a principio, esta deve ser a estrutura de vizinhança mais simples dentro das estruturas de diversificação

	float solution = 0.0;
	int sort1 = 0;
	int sort2 = 0;
	int aleatorio = 0;
	vector <int> mudarch;

	//sortear chaves para o VND
	while (sort1 == sort2)
	{
		aleatorio = rand() % (ac.numch_SIS - 1) + 1;
		sort1 = aleatorio; //primeiro sorteio

		aleatorio = rand() % (ac.numch_SIS - 1) + 1;
		sort2 = aleatorio; //segundo sorteio
	}

	mudarch.push_back(sort1);
	mudarch.push_back(sort2);

	solution = gvns.VND(mudarch);

	mudarch.clear();

	return solution;

}

float GVNS::v2_RVNS()
{
	//Descricao: escolher todas as chaves de um alimentador sorteado e aplicar o VND em todas as chaves alocadas nele

	int sort_AL = 0;
	int identificador = 0;
	vector <int> chavesAL;

	float solution = 0.0;

	// sorteando alimentador
	sort_AL = rand() % (num_AL - 1) + 1;

	//selecionando barras do alimentador
	identificador = 0;
	for (int i = 1; i < num_AL; i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			if (ac.chi[i][j] != 0 && ac.chf[i][j] != 0)
			{
				identificador++;

				if (sort_AL == i)
				{
					chavesAL.push_back(identificador);
				}
			}
		}
	}

	//solucao
	solution = gvns.VND(chavesAL);

	chavesAL.clear();

	return solution;
}

float GVNS::v3_RVNS()
{
	//Descricao: selecionar dois alimentadores e sortear uma chave para serem realocadas em posicoes aleatórias dentro do alimentador e fazer o VND em todas as chaves de cada alimentador

	float solution = 0.0;
	int sort_AL1 = 0;
	int sort_AL2 = 0;
	int sort_ch = 0;
	int numerochaves = 0;

	//copiar chaves e posicoes anteriores antes de usar aleatoriedade

	int chi_s[num_AL][linha_dados] = {};
	int chf_s[num_AL][linha_dados] = {};
	int ch_posicoes[num_AL][linha_dados] = {};

	for (int i = 1; i < num_AL; i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			chi_s[i][j] = ac.chi[i][j];
			chf_s[i][j] = ac.chf[i][j];
			ch_posicoes[i][j] = ac.posicaochaves[i][j];
		}
	}

	//outras variaveis

	int pos_aleat = 0;
	bool cond = false;
	bool repetir_passo = false;

	vector <int> posAL;
	vector <int> barsAL;
	vector <int> posCH;
	vector <int> mudarchaves;

	//sortear alimentadores
	sort_AL1 = 0;
	sort_AL2 = 0;
	while (sort_AL1 == sort_AL2)
	{
		sort_AL1 = rand() % (num_AL-1) + 1; //alimentador 1
		sort_AL2 = rand() % (num_AL-1) + 1; //alimentador 2
	}

	//para o alimentador sorteado

	repetir_passo = false;

novo_alimentador: //retorno para os alimentadores

	//identificar barras do alimentador
	for (int i = 1; i < linha_dados; i++)
	{
		if (ac.adjacente_chaves[sort_AL1][1][i] != 0)
		{
			barsAL.push_back(ac.adjacente_chaves[sort_AL1][1][i]); //pega as barras do alimentador
		}
	}

	//identificar posicoes possiveis de alocacao
	for (int i = 0; i < barsAL.size(); i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			if (ps.nof[j] == barsAL[i] && ps.candidato_aloc[j] == 1)
			{
				posAL.push_back(j); //pega a posicao das barras
			}
		}
	}

	//contar numero de chaves e guardar posicao delas
	numerochaves = 0;
	for (int i = 1; i < linha_dados; i++)
	{
		if (ac.chf[sort_AL1][i] != 0)
		{
			numerochaves++; //quantidade de chaves
			posCH.push_back(ac.posicaochaves[sort_AL1][i]); //guarda a posicao das chaves
		}
	}

	//sorteando chave
novosortv3:
	sort_ch = rand() % (numerochaves-1) + 1;

	//sorteando posicao aleatoria
	pos_aleat = rand() % (posAL.size()-1)+1;

	//analisando se ja tem esta chave
	cond = false;
	for (int i = 0; i < posCH.size(); i++)
	{
		if (posCH[i] == posAL[pos_aleat])
		{
			cond = true;
		}
	}

	//atribuindo chave a posicao aleatoria e selecionando todas do alimentador
	if (cond == false)
	{
		ac.posicaochaves[sort_AL1][sort_ch] = posAL[pos_aleat];
		ac.chi[sort_AL1][sort_ch] = ps.noi[posAL[pos_aleat]];
		ac.chf[sort_AL1][sort_ch] = ps.nof[posAL[pos_aleat]];
		ac.secoes_alimentador();

		//identificando chaves do alimentador
		numerochaves = 0;
		for (int i = 1; i < num_AL; i++)
		{
			for (int j = 1; j < linha_dados; j++)
			{
				if (ac.chi[i][j] != 0 && ac.chf[i][j] != 0)
				{
					numerochaves++;

					for (int k = 0; k < barsAL.size(); k++)
					{
						if (ac.chf[i][j] == barsAL[k])
						{
							mudarchaves.push_back(numerochaves);
						}
					}
				}
			}
		}

		//zerando vetores
		posAL.clear();
		barsAL.clear();
		posCH.clear();

		//retornar para repetir o processo para o segundo alimentador
		if (repetir_passo == false)
		{
			repetir_passo = true;
			sort_AL1 = sort_AL2;
			goto novo_alimentador;
		}
	}
	else
	{
		//ja tem esta chave no alimentador
		goto novosortv3;
	}

	//fazer VND
	solution = gvns.VND(mudarchaves);

	mudarchaves.clear();

	if (solution > gvns.incumbent_solution)
	{
		for (int i = 1; i < num_AL; i++)
		{
			for (int j = 1; j < linha_dados; j++)
			{
				ac.chi[i][j] = chi_s[i][j];
				ac.chf[i][j] = chf_s[i][j];
				ac.posicaochaves[i][j] = ch_posicoes[i][j];
			}
		}
	}

	return(solution);

}

float GVNS::v4_RVNS()
{
	//Descricao: seleciona todos os alimentadores e manter apenas uma chave de cada alimentador na posicao original, o restante deve ser realocado em posicoes aleatorias dentro do alimentador

	//copiar chaves e posicoes anteriores antes de usar aleatoriedade

	int chi_s[num_AL][linha_dados] = {};
	int chf_s[num_AL][linha_dados] = {};
	int ch_posicoes[num_AL][linha_dados] = {};

	for (int i = 1; i < num_AL; i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			chi_s[i][j] = ac.chi[i][j];
			chf_s[i][j] = ac.chf[i][j];
			ch_posicoes[i][j] = ac.posicaochaves[i][j];
		}
	}

	//outras variaveis 

	int savechi = 0;
	int savechf = 0;
	int saveposch = 0;

	int sortch = 0;
	bool dnv = false;

	int numberch = 0;
	vector <int> mudach;

	float solution = 0.0;

	for (int k = 1; k < num_AL; k++)
	{
		//sorteando chave para ser mantida
		sortch = rand() % (ac.numch_AL[k] - 1) + 1;

		saveposch = ac.posicaochaves[k][sortch];
		savechi = ac.chi[k][sortch];
		savechf = ac.chf[k][sortch];

		ps.leitura_parametros();
		fxp.fluxo_potencia();
		sorteiochaves(ac.numch_AL[k], fxp.camadaAL[k], ac.posicaochaves[k], alimentadores[k]);

		for (int i = 1; i < num_AL; i++)
		{
			for (int j = 1; j < linha_dados; j++)
			{
				ac.chi[i][j] = ps.noi[ac.posicaochaves[i][j]];
				ac.chf[i][j] = ps.nof[ac.posicaochaves[i][j]];
			}
		}

		//conservando a chave
		ac.posicaochaves[k][sortch] = saveposch;
		ac.chi[k][sortch] = savechi;
		ac.chf[k][sortch] = savechf;
		ac.secoes_alimentador();

		//salvando as chaves para o vnd
		numberch = 0;
		for (int i = 1; i < num_AL; i++)
		{
			for (int j = 1; j < linha_dados; j++)
			{
				if (ac.chf[i][j] != 0 && ac.chi[i][j] != 0)
				{
					numberch++;
					if (i == k)
					{
						mudach.push_back(numberch);
					}
				}
			}
		}
	}

	solution = gvns.VND(mudach);

	if (solution > gvns.incumbent_solution)
	{
		for (int i = 1; i < num_AL; i++)
		{
			for (int j = 1; j < linha_dados; j++)
			{
				ac.chi[i][j] = chi_s[i][j];
				ac.chf[i][j] = chf_s[i][j];
				ac.posicaochaves[i][j] = ch_posicoes[i][j];
			}
		}
	}

	mudach.clear();

	return solution;

}

//############################################################################################

int main()
{
	srand(static_cast <unsigned int> (time(NULL)));	//faz a aleatoriedade com base no relogio

	int itGVNS = 0;

	//variaveis a serem analisadas
	gvns.current_solution = 0.0;
	gvns.incumbent_solution = 0.0;

	gvns.q_rvns1 = 0;
	gvns.q_rvns2 = 0;
	gvns.q_rvns3 = 0;
	gvns.q_rvns4 = 0;
	gvns.q_vnd1 = 0;
	gvns.q_vnd2 = 0;

	//faz a leitura dos parametros do circuito
	ps.leitura_parametros();

	//somatorio da potencia total do sistema
	ps.somatorio_potencia();

	//resolve o fluxo de potencia
	fxp.fluxo_potencia();

	//tirando de pu, para conferir - imprime o fluxo de potencia
	//fxp.valores_nominais_tensao();

	//define quantas chaves serao alocadas em cada alimentador
	ac.criterio_numero_de_chaves();

	//primeira alocacao: esta eh feita de forma aleatoria
	gvns.primeiraaloc();

	//iniciando do GVNS

	//primeiras chaves
	for (int i = 1; i < num_AL; i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			ac.chi[i][j] = ps.noi[ac.posicaochaves[i][j]];
			ac.chf[i][j] = ps.nof[ac.posicaochaves[i][j]];
		}
	}

	//contado total de chaves
	for (int i = 1; i < num_AL; i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			if (ac.chi[i][j] != 0 && ac.chf[i][j] != 0)
			{
				ac.numch_SIS++;
			}
		}
	}

	//imprimindo chaves iniciais:
	cout << "Chaves Iniciais:" << endl;

	for (int i = 1; i < num_AL; i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			if (ac.chi[i][j] != 0 && ac.chf[i][j] != 0)
			{
				cout << ac.chi[i][j];
				cout << "|--|";
				cout << ac.chf[i][j] << endl;
			}
		}
	}
	cout << "\n";

	ps.leitura_parametros();
	fxp.fluxo_potencia();
	ac.secoes_alimentador();

	for (int i = 1; i < num_AL; i++)
	{
		gvns.analisepos = i;
		gvns.incumbent_solution = ac.calculo_funcao_objetivo();
	}
	
	gvns.current_solution = gvns.incumbent_solution;

	cout << "$ Sol. Inicial: " << gvns.incumbent_solution << endl;

nmgvns:

	itGVNS++;

metaheuristicGVNS:

	cout << "# 1" << endl;
	gvns.current_solution = gvns.v1_RVNS();

	if (gvns.current_solution < gvns.incumbent_solution)
	{
		gvns.incumbent_solution = gvns.current_solution;
		gvns.q_rvns1++;

		cout << "# 1-RVNS: " << gvns.incumbent_solution << endl;

		goto metaheuristicGVNS;
	}

	cout << "# 2" << endl;
	gvns.current_solution = gvns.v2_RVNS();

	if (gvns.current_solution < gvns.incumbent_solution)
	{
		gvns.incumbent_solution = gvns.current_solution;
		gvns.q_rvns2++;

		cout << "# 2-RVNS: " << gvns.incumbent_solution << endl;
		
		goto metaheuristicGVNS;
	}

	cout << "# 3" << endl;
	gvns.current_solution = gvns.v3_RVNS();

	if (gvns.current_solution < gvns.incumbent_solution)
	{
		gvns.incumbent_solution = gvns.current_solution;
		gvns.q_rvns3++;

		cout << "# 3-RVNS: " << gvns.incumbent_solution << endl;

		goto metaheuristicGVNS;
	}

	cout << "# 4" << endl;
	gvns.current_solution = gvns.v4_RVNS();

	if (gvns.current_solution < gvns.incumbent_solution)
	{
		gvns.incumbent_solution = gvns.current_solution;
		gvns.q_rvns4++;

		cout << "# 4-RVNS: " << gvns.incumbent_solution << endl;

		goto metaheuristicGVNS;
	}

	if(itGVNS < 10) { goto nmgvns; }

	////////////////////////////////////////////////////////////
	//Fim GVNS

	cout << "----------------------------------" << endl;

	//chaves finais
	cout << "Chaves Finais:" << endl;

	for (int i = 1; i < num_AL; i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			if (ac.chi[i][j] != 0 && ac.chf[i][j] != 0)
			{
				cout << ac.chi[i][j];
				cout << " |--| ";
				cout << ac.chf[i][j] << endl;
			}
		}
	}
	cout << "\n";

	//total iteracoes
	cout << "Numero de iteracoes: " << itGVNS << endl;
	cout << "\n";

	cout << "Operacoes nas vizinhacas:" << endl;
	cout << "rvns 1: " << gvns.q_rvns1 << endl;
	cout << "rvns 2: " << gvns.q_rvns2 << endl;
	cout << "rvns 3: " << gvns.q_rvns3 << endl;
	cout << "rvns 4: " << gvns.q_rvns4 << endl;
	cout << "vnd 1: " << gvns.q_vnd1 << endl;
	cout << "vnd 2: " << gvns.q_vnd2 << endl;
	cout << "\n";

	//fluxo de potencia:
	cout << "Numero de fluxo de potencia da simulacao: " << fxp.contadorFXP << endl;
	cout << "\n";

	//Solucao final
	cout << "Solucao otima: " << gvns.incumbent_solution << endl;
	cout << "\n\n";

	cout << "------------------------------------------------------------------------------------------------------------------" << endl;
	cout << "Fim" << endl;

	return 0;
}