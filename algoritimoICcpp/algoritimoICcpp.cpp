//IC 2020
//ALOCAÇÃO DE CHAVES DE MANOBRA PARA RESTAURAÇÃO

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <cmath>
#include <complex>
#include <vector>
#include <algorithm>
#include <tuple>

using namespace std;

// DADOS DO SISTEMA DE 136 BARRAS ---------------------------------------

//Dados para o arquivo txt
#define linha_dados 157 //numero de linhas da matriz de dados +1
#define coluna_dados 9 //numero de colunas da matriz de dados +1

//Caracteristicas Alimentadores e Subestacoes ----------------------------

#define num_AL 9 //numero de alimentadores +1 por conta do cpp
#define estado_restaurativo_pu 0.93 //tensao minima no estado restaurativo
#define capSE 100000000 //potencia maxima de cada alimentador

int alimentadores[num_AL] = { 0, 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007 }; // alimentadores

//Chave a cada quantos kVA?
#define parametroCH_kVA 800;

//PARAMETROS FUNCAO OBJETIVO E METAHEURISTICA ------------------------------------

#define tempo_falha 4 //numero de horas que o sistema fica em estado restaurativo
#define tempo_isolacao 0.12 //tempo necessario para fazer as manobras em horas
#define taxa_falhas 0.18 //taxa de falhas por km no ano
#define custoKWh 0.12 // em real 0.53187 (cotação 2017 ANEEL - Elektro - Sudeste)
#define criterio_parada 3 //criterio de parada da metaheuristica
#define numero_simulacoes 26 //numero de simulacoes que o algoritmo faz

//Caracteristicas Fluxo de Potencia ------------------------------------

//valores base, ou, de referencia
float sref = 100 * pow(10, 3); //100MVA
float vref = 13800; //13.8kV
float zref = (vref * vref) / sref;
float iref = sref / vref;

//tensao inicial nos nós do sistema: vref * (valor para dar a tensao inicial)
float pu_inicial_alimentador = 1.05;
float tensao_inicial_nos = (vref)*pu_inicial_alimentador;

//critério de convergencia do fluxo de potencia
complex <float> criterio_conv = 1 * pow(10, -4);
float epsilon = abs(criterio_conv);

int max_interacao = 8;


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
	int estado_swt_vanila[linha_dados] = {}; //estado das chaves para o sistema vanila

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

private:

	void camadas(int alimentador, int camadaalimentador[linha_dados][linha_dados]);
	void backward_sweep(int alimentador, int almt);
	void forward_sweep(int alimentador, int alt);

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
	int antpos[num_AL][linha_dados] = {}; //posicao da chave no sistema

	float fo_al[num_AL] = {};
	float fo_al_save[num_AL] = {};

	void chaves_anteriores();
	void volta_chaves_anteriores();

	void criterio_numero_de_chaves();
	void secoes_alimentador();
	float calculo_funcao_objetivo(int p_AL);
	float calculo_funcao_objetivo_geral();

private:

	tuple <int, float> contagem_criterio(int camada[linha_dados][linha_dados]); //criterio para a contagem de quantas chaves alocar em cada alimentador do sistema teste
	void adjacentes(int posicao[linha_dados], int adj[linha_dados][linha_dados], int alimentador); //calcula os adjacentes das chaves e da secao do alimentador
	float FO(float potencia_secao, float comprimeto_secao, float ens_isolacao);
	float energia_suprida(int AL, vector<int>barras_AL);

}ac;

class GVNS
{
public:

	int antchivnd[num_AL][linha_dados] = {}; //barra inicial da chave
	int antchfvnd[num_AL][linha_dados] = {}; //barra final da chave
	int antposvnd[num_AL][linha_dados] = {}; //posicao da chave no sistema

	float fo_al_savevnd[num_AL] = {};

	void chaves_anterioresVND(int alm);
	void volta_chaves_anterioresVND(int almv);
	void fo_anteriorVND(int almfo);
	void volta_fo_anteriorVND(int almfov);

	void primeiraaloc();
	void sorteiochaves(int numch, int camada[linha_dados][linha_dados], int posicao_camada[linha_dados], int alimentador); //sorteio inicial das chaves

}gvns;

class RVNS : public GVNS
{
public:

	int q_rvns1 = 0;
	int q_rvns2 = 0;
	int q_rvns3 = 0;
	int q_rvns4 = 0;

	float v1_RVNS(float incumbentmainv1);
	float v2_RVNS(float incumbentmainv2);
	float v3_RVNS(float incumbentmainv3);
	float v4_RVNS(float incumbentmainv4);
}rvns;

class VND : public GVNS
{
public:

	int q_vnd1 = 0;
	int q_vnd2 = 0;

	float VND_intensificacao(int ch, float sol_incumbent);
	float v1_VND(int ch1, float incumbentv1); //mover para adjacente
	float v2_VND(int ch2, float incumbentv2); //mover para adjacente do adjacente

}vnd;

//------------------------------------------------------------

void ParametrosSistema::leitura_parametros()
{
	FILE* arquivo;

	if ((arquivo = fopen("dados136.txt", "r")) == NULL)
		return;

	for (int i = 1; i < linha_dados; i++)
	{
		fscanf(arquivo, "%d%d%f%f%f%f%d%d%f", &ps.noi[i], &ps.nof[i], &ps.lt_r[i], &ps.lt_x[i], &ps.s_nofr[i], &ps.s_nofq[i], &ps.candidato_aloc[i], &ps.estado_swt_vanila[i], &ps.dist_no[i]);
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


	//estado das chaves
	for (int k = 1; k < linha_dados; k++)
	{
		ps.estado_swt[k] = ps.estado_swt_vanila[k];
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

void FluxoPotencia::backward_sweep(int alimentador, int almt)
{
	int posi = 0;
	int posf = 0;

	//atribuição das correntes nas barras
	for (int i = 1; i < linha_dados; i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			for (int k = 1; k < linha_dados; k++)
			{
				if (fxp.camadaAL[almt][i][j] == ps.nof[k])
				{
					corrente_pu[k] = conj(ps.pu_s_nof[k] / tensao_pu[k]);
				}
			}
		}
	}

	// somatorio das correntes nos ramos
	for (int i = linha_dados - 1; i > 2; i--) //o i aqui vai ate 1 pq a linha 1 é a do alimentador e i-1 sera a linha 1
	{
		for (int k = linha_dados - 1; k > 0; k--)
		{
			if (fxp.camadaAL[almt][i][k] != 0)
			{
				for (int r = 1; r < linha_dados; r++)
				{
					if (fxp.camadaAL[almt][i - 1][r] == 0) { continue; }

					for (int j = 1; j < linha_dados; j++)
					{
						if (fxp.camadaAL[almt][i][k] == ps.nof[j] && fxp.camadaAL[almt][i - 1][r] == ps.noi[j])
						{
							posf = 0;
							posi = 0;

							for (int o = 1; o < linha_dados; o++)
							{
								if (ps.nof[o] == fxp.camadaAL[almt][i - 1][r] && ps.estado_swt[o] == 1)
								{
									if (ps.candidato_aloc[o] == 1 || ps.noi[o] == alimentador)
									{
										posf = o;
										break;
									}
								}
							}

							for (int o = 1; o < linha_dados; o++)
							{
								if (ps.nof[o] == fxp.camadaAL[almt][i][k] && ps.estado_swt[o] == 1)
								{
									if (ps.candidato_aloc[o] == 1 || ps.noi[o] == alimentador)
									{
										posi = o;
										break;
									}
								}
							}

							fxp.corrente_pu[posf] = fxp.corrente_pu[posf] + fxp.corrente_pu[posi];
						}
						else if (fxp.camadaAL[almt][i][k] == ps.noi[j] && fxp.camadaAL[almt][i - 1][r] == ps.nof[j])
						{
							posf = 0;
							posi = 0;

							for (int o = 1; o < linha_dados; o++)
							{
								if (ps.nof[o] == fxp.camadaAL[almt][i - 1][r] && ps.estado_swt[o] == 1)
								{
									if (ps.candidato_aloc[o] == 1 || ps.noi[o] == alimentador)
									{
										posf = o;
										break;
									}
								}
							}

							for (int o = 1; o < linha_dados; o++)
							{
								if (ps.nof[o] == fxp.camadaAL[almt][i][k] && ps.estado_swt[o] == 1)
								{
									if (ps.candidato_aloc[o] == 1 || ps.noi[o] == alimentador)
									{
										posi = o;
										break;
									}
								}
							}

							fxp.corrente_pu[posf] = fxp.corrente_pu[posf] + fxp.corrente_pu[posi];
						}
					}
				}
				
			}
		}
	}
	
}

void FluxoPotencia::forward_sweep(int alimentador, int alt)
{

	complex <float> unit = fxp.tensao_inicial; //complexo unitario em pu para a tensao
	int vf = 0;
	int vi = 0;

	//atribuindo a tensao para o restante dos nós

	for (int i = 1; i < (linha_dados - 2); i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			if (fxp.camadaAL[alt][i][j] == 0) { continue; } //continua se for 0 para adiantar o processo

			for (int k = 1; k < linha_dados; k++)
			{
				if (fxp.camadaAL[alt][i + 1][k] == 0) { continue; }

				for (int t = 1; t < linha_dados; t++)
				{
					if (fxp.camadaAL[alt][i][j] == alimentador)
					{
						vf = 0;
						for (int r = 1; r < linha_dados; r++)
						{
							if (ps.noi[r] == alimentador)
							{
								vf = r;
								break;
							}
						}

						fxp.tensao_pu[vf] = unit - (fxp.corrente_pu[vf] * ps.pu_lt[vf]);
						break;
					}
					else if (fxp.camadaAL[alt][i][j] == ps.noi[t] && fxp.camadaAL[alt][i + 1][k] == ps.nof[t] && fxp.camadaAL[alt][i][j] != alimentador)
					{
						vf = 0;
						vi = 0;

						for (int r = 1; r < linha_dados; r++)
						{
							if (ps.nof[r] == fxp.camadaAL[alt][i][j] && ps.estado_swt[r] == 1)
							{
								vi = r;
								break;
							}
						}

						for (int r = 1; r < linha_dados; r++)
						{
							if (ps.nof[r] == fxp.camadaAL[alt][i + 1][k] && ps.estado_swt[r] == 1)
							{
								vf = r;
								break;
							}
						}

						fxp.tensao_pu[vf] = fxp.tensao_pu[vi] - (fxp.corrente_pu[vf] * ps.pu_lt[t]);

					}

					else if (fxp.camadaAL[alt][i][j] == ps.nof[t] && fxp.camadaAL[alt][i + 1][k] == ps.noi[t] && fxp.camadaAL[alt][i][j] != alimentador)
					{
						vf = 0;
						vi = 0;

						for (int r = 1; r < linha_dados; r++)
						{
							if (ps.nof[r] == fxp.camadaAL[alt][i][j] && ps.estado_swt[r] == 1)
							{
								vi = r;
								break;
							}
						}

						for (int r = 1; r < linha_dados; r++)
						{
							if (ps.nof[r] == fxp.camadaAL[alt][i + 1][k] && ps.estado_swt[r] == 1)
							{
								vf = r;
								break;
							}
						}

						fxp.tensao_pu[vf] = fxp.tensao_pu[vi] - (fxp.corrente_pu[vf] * ps.pu_lt[t]);
					}
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
		cout << ((arg(fxp.tensao_pu[i])) * 180.0 / 3.14) << "\n";

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

	complex <float> convergencia_tensao[linha_dados];

	int som = 0;
	int iteracao = 0;

	//valores iniciais
	for (int i = 1; i < linha_dados; i++)
	{
		tensao_aux[i] = tensao_inicial;
		corrente_aux[i] = 0.0;

		fxp.tensao_pu[i] = tensao_inicial;
		fxp.corrente_pu[i] = 0.0;
	}
	//definir as camadas de fxp.camadas de cada alimentador

	for (int i = 1; i < num_AL; i++)
	{
		camadas(alimentadores[i], fxp.camadaAL[i]);
	}

	//estado das chaves
	for (int k = 1; k < linha_dados; k++)
	{
		ps.estado_swt[k] = ps.estado_swt_vanila[k];
	}

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
		}

		//1 passo: BACKWARD

		//corrente nos ramos
		for (int i = 1; i < num_AL; i++)
		{
			backward_sweep(alimentadores[i], i);
		}

		//2 passo: FORWARD

		//tensoes nodais
		for (int i = 1; i < num_AL; i++)
		{
			forward_sweep(alimentadores[i], i);
		}

		//comparação de critério satisfeito

		for (int k = 1; k < linha_dados; k++)
		{
			convergencia_tensao[k] = fxp.tensao_pu[k] - tensao_aux[k];
		}

		for (int k = 1; k < linha_dados; k++)
		{
			if (abs(convergencia_tensao[k]) >= epsilon)
			{
				som += 1;
			}
		}

		if (som == 0 /*|| iteracao == max_interacao*/)
		{
			criterio_satisfeito = true;
		}

	}

	//cout << " max it:     " << iteracao << endl;

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

	num_crit = floor(num);
	//num_crit = round(num);

	if (num_crit < 2) { num_crit = 2; }

	//num_crit = 2;

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

float AlocacaoChaves::energia_suprida(int AL, vector<int>barras_AL)
{
	bool analise = false;
	bool conta_pot = false;
	float energia_sup = 0.0;
	complex <float> capacidadeSE;

	// ANALISANDO FLUXO DE POTENCIA:

	//faz o fluxo de potencia para a situação do remanejamento
	fxp.fluxo_potencia();

	analise = true; //verdadeiro se atende as condiçoes, falso caso contrario - assume-se inicialmente

	// 1. Níveis de tensao e corrente
	for (int y = 1; y < linha_dados; y++)
	{
		//para a tensao:
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

		if (abs(capacidadeSE) > capSE)
		{
			analise = false;
			break;
		}
	}


	//feita a analise do fluxo de putencia:
	if (analise == true)
	{
		energia_sup = 0.0;

		//ver a energia suprida neste caso
		for (int i = 0; i < barras_AL.size(); i++)
		{
			conta_pot = false;

			for (int u = 1; u < num_AL; u++)
			{
				if (AL == u) { continue; }
				else
				{
					for (int k = 2; k < linha_dados; k++) //começa no 2 para tirar o al da analise
					{
						for (int j = 1; j < linha_dados; j++)
						{
							if (fxp.camadaAL[u][k][j] != 0)
							{
								if (fxp.camadaAL[u][k][j] == barras_AL[i]) { conta_pot = true; }
							}
						}
					}
				}
			}

			if (conta_pot == true)
			{
				for (int k = 1; k < linha_dados; k++)
				{
					if (barras_AL[i] == ps.nof[k])
					{
						energia_sup = energia_sup + ps.s_nofr[k];
					}
				}
			}
		}
	}
	else
	{
		energia_sup = 0.0;
	}

	
	//retornar valor a variavel
	return energia_sup;
}

float AlocacaoChaves::FO(float potencia_secao, float comprimento, float ens_isolacao)
{
	float resultado = 0.0;

	resultado = 0.0;

	resultado = (ens_isolacao * tempo_isolacao * custoKWh);
	resultado = (potencia_secao * custoKWh * tempo_falha) + resultado;
	resultado = taxa_falhas * comprimento * resultado;

	return resultado;
}

void AlocacaoChaves::chaves_anteriores()
{
	//salva as chaves
	for (int i = 1; i < num_AL; i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			ac.antchi[i][j] = ac.chi[i][j];
			ac.antchf[i][j] = ac.chf[i][j];
			ac.antpos[i][j] = ac.posicaochaves[i][j];
		}
	}
}

void AlocacaoChaves::volta_chaves_anteriores()
{
	//volta as chaves
	for (int i = 1; i < num_AL; i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			ac.chi[i][j] = ac.antchi[i][j];
			ac.chf[i][j] = ac.antchf[i][j];
			ac.posicaochaves[i][j] = ac.antpos[i][j];
		}
	}
}

float AlocacaoChaves::calculo_funcao_objetivo(int p_AL)
{
	float comprimento_secao = 0.0;
	float potencia_W = 0.0;
	float valorFO = 0.0;
	float chamadaFO = 0.0;
	float potencia_isolacao = 0.0;
	float ENSotima = 0.0;
	float ens = 0.0;
	float resultado_FO = 0.0;

	int contcondicao = 0;

	vector <int> secao = {};
	vector <int> posicao = {};
	vector <int> barras = {};
	vector <vector<int>> analise_remanejamento = {};
	vector <vector<int>> remanej_cargas = {};

	posicao.clear();
	secao.clear();
	analise_remanejamento.clear();
	remanej_cargas.clear();

	//deve-se analisar todas as secoes para os valores da funcao obj

	ac.secoes_alimentador();

	valorFO = 0.0;

	ps.leitura_parametros();
	fxp.fluxo_potencia();
	
	//barras do alimentador analizado
	barras.clear();
	for (int i = 1; i < linha_dados; i++)
	{
		if (ac.adjacente_chaves[p_AL][1][i] != 0)
		{
			barras.push_back(ac.adjacente_chaves[p_AL][1][i]);
		}
	}

	//secao j do alimentador p_AL
	for (int j = 1; j < linha_dados; j++)
	{
		//ver se vale a pena fazer o laço, se o vetor estiver zerado é só custo computacional a toa
		contcondicao = 0;

		for (int k = 0; k < linha_dados; k++)
		{
			if (ac.secoes_chaves[p_AL][j][k] != 0)
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
				if (ac.secoes_chaves[p_AL][j][k] == ps.nof[y])
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
				if (ac.adjacente_chaves[p_AL][1][k] == ps.nof[y])
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
				if (ac.secoes_chaves[p_AL][j][k] == 0 || ac.posicaochaves[p_AL][y] == 0) { continue; }
				else if (ac.secoes_chaves[p_AL][j][k] == ps.nof[ac.posicaochaves[p_AL][y]] || ac.secoes_chaves[p_AL][j][k] == ps.noi[ac.posicaochaves[p_AL][y]])
				{
					ps.estado_swt[ac.posicaochaves[p_AL][y]] = 0;
				}
			}
		}

		if (j == 1) //desligar o disjuntor da subestacao
		{
			for (int k = 1; k < linha_dados; k++)
			{
				if (ps.noi[k] == alimentadores[p_AL])
				{
					ps.estado_swt[k] = 0;
				}
			}
		}
		
		// 2b) cenario da falta
		//zerar falha na camada
		for (int k = 1; k < linha_dados; k++)
		{
			for (int y = 1; y < linha_dados; y++)
			{
				for (int t = 1; t < linha_dados; t++)
				{
					if (fxp.camadaAL[p_AL][k][y] == ac.secoes_chaves[p_AL][j][t])
					{
						fxp.camadaAL[p_AL][k][y] = 0;
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
				if (fxp.camadaAL[p_AL][k][y] != 0)
				{
					secao.push_back(fxp.camadaAL[p_AL][k][y]);
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
							if (secao[t] == ac.secoes_chaves[p_AL][m][n])
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
				if (ac.adjacente_chaves[p_AL][posicao[k]][t] != 0)
				{
					secao.push_back(ac.adjacente_chaves[p_AL][posicao[k]][t]);
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
			ens = ps.potencia_al[p_AL];

			//somente uma opcao para remanejamento de cada vez
			for (int a = 0; a < remanej_cargas.size(); a++)
			{
				potencia_W = 0.0;

				for (int b = 0; b < remanej_cargas[a].size(); b++)
				{
					//fechando chave da secao
					for (int k = 1; k < linha_dados; k++)
					{
						for (int y = 1; y < linha_dados; y++)
						{
							if (ac.secoes_chaves[p_AL][j][k] == 0 || ac.posicaochaves[p_AL][y] == 0) { continue; }
							else if (ac.secoes_chaves[p_AL][j][k] == ps.nof[ac.posicaochaves[p_AL][y]] || ac.secoes_chaves[p_AL][j][k] == ps.noi[ac.posicaochaves[p_AL][y]])
							{
								ps.estado_swt[ac.posicaochaves[p_AL][y]] = 0;
							}
						}
					}

					//abrindo chave
					ps.estado_swt[remanej_cargas[a][b]] = 1;

					potencia_W = ac.energia_suprida(p_AL, barras);

					ps.estado_swt[remanej_cargas[a][b]] = 0;

					if (potencia_W != 0) { break; }

				}

				ens = ens - potencia_W;
			}

			if (ens == ps.potencia_al[p_AL])
			{
				//nao tem como fazer manobra, a ENS será os adjacentes da chave
				ens = 0.0;

				for (int y = 1; y < linha_dados; y++)
				{
					for (int h = 1; h < linha_dados; h++)
					{
						if (ps.nof[h] == ac.adjacente_chaves[p_AL][j][y])
						{
							ens = ens + ps.s_nofr[h];
						}
					}
				}
			}

			ps.leitura_parametros();
			fxp.fluxo_potencia();
		}
		else
		{
			//nao tem como fazer manobra, a ENS será os adjacentes da chave
			ens = 0.0;

			for (int y = 1; y < linha_dados; y++)
			{
				for (int h = 1; h < linha_dados; h++)
				{
					if (ps.nof[h] == ac.adjacente_chaves[p_AL][j][y])
					{
						ens = ens + ps.s_nofr[h];
					}
				}
			}

			ps.leitura_parametros();
			fxp.fluxo_potencia();
		}

		analise_remanejamento.clear();
		remanej_cargas.clear();
		secao.clear();
		posicao.clear();


		// 4) chamando a funcao objetivo
		chamadaFO = 0.0;

		chamadaFO = FO(ens, comprimento_secao, potencia_isolacao);

		valorFO += chamadaFO;
	}

	//zerar vetor barras
	barras.clear();

	//analisando melhor caso
	if (valorFO < fo_al[p_AL])
	{
		fo_al[p_AL] = valorFO;
	}
	else
	{
		volta_chaves_anteriores();
		ac.secoes_alimentador();
	}

	//calculo F0 geral
	resultado_FO = 0.0;
	for (int i = 1; i < num_AL; i++)
	{
		resultado_FO += fo_al[i];
	}

	return resultado_FO;
}

float AlocacaoChaves::calculo_funcao_objetivo_geral()
{
	float comprimento_secao = 0.0;
	float potencia_W = 0.0;
	float valorFO = 0.0;
	float chamadaFO = 0.0;
	float potencia_isolacao = 0.0;
	float ENSotima = 0.0;
	float ens = 0.0;
	float ens2 = 0.0;
	float resultado_FO = 0.0;

	int contcondicao = 0;

	vector <int> secao = {};
	vector <int> posicao = {};
	vector <int> barras = {};
	vector <vector<int>> analise_remanejamento = {};
	vector <vector<int>> remanej_cargas = {};
	vector <vector<int>> remanej_cargas2 = {};

	posicao.clear();
	secao.clear();
	analise_remanejamento.clear();
	remanej_cargas.clear();

	//deve-se analisar todas as secoes para os valores da funcao obj

	ac.secoes_alimentador();

	valorFO = 0.0;

	ps.leitura_parametros();
	fxp.fluxo_potencia();

	//secao j do alimentador w

	for (int w = 1; w < num_AL; w++)
	{
		//barras do alimentador w
		barras.clear();
		for (int q = 1; q < linha_dados; q++)
		{
			if (ac.adjacente_chaves[w][1][q] != 0)
			{
				barras.push_back(ac.adjacente_chaves[w][1][q]);
			}
		}

		//secao j do alimetador w
		for (int j = 1; j < linha_dados; j++)
		{
			//ver se vale a pena fazer o laço, se o vetor estiver zerado é só custo computacional a toa
			contcondicao = 0;

			for (int k = 0; k < linha_dados; k++)
			{
				if (ac.secoes_chaves[w][j][k] != 0)
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
					if (ac.secoes_chaves[w][j][k] == ps.nof[y])
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
					if (ac.adjacente_chaves[w][1][k] == ps.nof[y])
					{
						potencia_isolacao = potencia_isolacao + ps.s_nofr[y];
					}
				}
			}


			// 2) agora deve-se fazer o devido chaveamento
			// 2a) isolando secao j abrindo a chave
			for (int k = 1; k < linha_dados; k++)
			{
				for (int y = 1; y < linha_dados; y++)
				{
					if (ac.secoes_chaves[w][j][k] == 0 || ac.posicaochaves[w][y] == 0) { continue; }
					else if (ac.secoes_chaves[w][j][k] == ps.nof[ac.posicaochaves[w][y]] || ac.secoes_chaves[w][j][k] == ps.noi[ac.posicaochaves[w][y]])
					{
						ps.estado_swt[ac.posicaochaves[w][y]] = 0;
					}
				}
			}

			if (j == 1) //desligar o disjuntor da subestacao
			{
				for (int k = 1; k < linha_dados; k++)
				{
					if (ps.noi[k] == alimentadores[w])
					{
						ps.estado_swt[k] = 0;
					}
				}
			}

			// 2b) cenario da falta
			//zerar falha na camada
			for (int k = 1; k < linha_dados; k++)
			{
				for (int y = 1; y < linha_dados; y++)
				{
					for (int t = 1; t < linha_dados; t++)
					{
						if (fxp.camadaAL[w][k][y] == ac.secoes_chaves[w][j][t])
						{
							fxp.camadaAL[w][k][y] = 0;
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
					if (fxp.camadaAL[w][k][y] != 0)
					{
						secao.push_back(fxp.camadaAL[w][k][y]);
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
								if (secao[t] == ac.secoes_chaves[w][m][n])
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
					if (ac.adjacente_chaves[w][posicao[k]][t] != 0)
					{
						secao.push_back(ac.adjacente_chaves[w][posicao[k]][t]);
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
				ens = ps.potencia_al[w];

				//somente uma opcao para remanejamento de cada vez, no caso 1, uma chave para remanejamento é aberta de cada vez
				for (int a = 0; a < remanej_cargas.size(); a++)
				{
					potencia_W = 0.0;

					for (int b = 0; b < remanej_cargas[a].size(); b++)
					{
						//abrindo chaves
						for (int k = 1; k < linha_dados; k++)
						{
							for (int y = 1; y < linha_dados; y++)
							{
								if (ac.secoes_chaves[w][j][k] == 0 || ac.posicaochaves[w][y] == 0) { continue; }
								else if (ac.secoes_chaves[w][j][k] == ps.nof[ac.posicaochaves[w][y]] || ac.secoes_chaves[w][j][k] == ps.noi[ac.posicaochaves[w][y]])
								{
									ps.estado_swt[ac.posicaochaves[w][y]] = 0;
								}
							}
						}

						ps.estado_swt[remanej_cargas[a][b]] = 1; //fecha chave remanejamento

						potencia_W = ac.energia_suprida(w, barras);

						ps.estado_swt[remanej_cargas[a][b]] = 0; //abre chave remanejamento

						if (potencia_W != 0)
						{
							break;
						}

					}

					ens = ens - potencia_W;
				}

				if (ens == ps.potencia_al[w])
				{
					//nao tem como fazer manobra, a ENS será os adjacentes da chave
					ens = 0.0;

					for (int y = 1; y < linha_dados; y++)
					{
						for (int h = 1; h < linha_dados; h++)
						{
							if (ps.nof[h] == ac.adjacente_chaves[w][j][y])
							{
								ens = ens + ps.s_nofr[h];
							}
						}
					}
				}

				ps.leitura_parametros();
				fxp.fluxo_potencia();
			}
			else
			{
				//nao tem como fazer manobra, a ENS será os adjacentes da chave
				ens = 0.0;

				for (int y = 1; y < linha_dados; y++)
				{
					for (int h = 1; h < linha_dados; h++)
					{
						if (ps.nof[h] == ac.adjacente_chaves[w][j][y])
						{
							ens = ens + ps.s_nofr[h];
						}
					}
				}

				ps.leitura_parametros();
				fxp.fluxo_potencia();
			}

			analise_remanejamento.clear();
			remanej_cargas.clear();
			secao.clear();
			posicao.clear();


			// 4) chamando a funcao objetivo
			chamadaFO = 0.0;

			chamadaFO = FO(ens, comprimento_secao, potencia_isolacao);

			valorFO += chamadaFO;
		}

		fo_al[w] = valorFO;
		valorFO = 0.0;
		barras.clear();
	}
	
	valorFO = 0.0;
	for (int i = 1; i < num_AL; i++)
	{
		valorFO = valorFO + fo_al[i];
	}

	return valorFO;
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

void GVNS::chaves_anterioresVND(int alm)
{
	//salva as chaves do alimentador informado
	for (int j = 1; j < linha_dados; j++)
	{
		antchivnd[alm][j] = ac.chi[alm][j];
		antchfvnd[alm][j] = ac.chf[alm][j];
		antposvnd[alm][j] = ac.posicaochaves[alm][j];
	}
}

void GVNS::volta_chaves_anterioresVND(int almv)
{
	//volta as chaves	 
	for (int j = 1; j < linha_dados; j++)
	{
		ac.chi[almv][j] = antchivnd[almv][j];
		ac.chf[almv][j] = antchfvnd[almv][j];
		ac.posicaochaves[almv][j] = antposvnd[almv][j];
	}
}

void GVNS::fo_anteriorVND(int almfo)
{
	//salva os valores da funcao objetivo anterior
	fo_al_savevnd[almfo] = ac.fo_al[almfo];
}

void GVNS::volta_fo_anteriorVND(int almfov)
{
	//volta os valores para o vetor da funcao objetivo
	ac.fo_al[almfov] = fo_al_savevnd[almfov];
}

float VND::v1_VND(int ch1, float incumbentv1)
{
	//troca as chaves para o vizinho, ou seja, para o adjacente
	//chavesv1 é o numero das chaves

	int identf = 0;
	vector <int> pos_vizinhas;

	bool repetepos = false;
	int pos_vetorcaso_nao_melhor = 0;

	float soluc = 0.0;
	float mainsoluc = 0.0;

	vector<int>auxfunc; //vetor auxiliar das funcoes
	vector <vector<int>> identch; //sempre tera duas posicoes:: 1)chi - 2)chf

	mainsoluc = incumbentv1;

	// 1) localizar chaves sorteadas

	identf = 0;

	for (int j = 1; j < num_AL; j++)
	{
		for (int k = 1; k < linha_dados; k++)
		{
			if (ac.chi[j][k] != 0 && ac.chf[j][k] != 0)
			{
				identf++;

				if (identf == ch1)
				{
					auxfunc.push_back(ac.chi[j][k]);
					auxfunc.push_back(ac.chf[j][k]);

					identch.push_back(auxfunc);
					auxfunc.clear();
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
		soluc = incumbentv1;
		pos_vetorcaso_nao_melhor = ac.posicaochaves[pos1ch][pos2ch];
		bool igual_pos = false;

		for (int e = 0; e < pos_vizinhas.size(); e++)
		{
			igual_pos = false;
			for (int d = 1; d < linha_dados; d++)
			{
				if (ac.posicaochaves[pos1ch][d] != 0 && ac.posicaochaves[pos1ch][d] == pos_vizinhas[e])
				{
					igual_pos = true;
				}
			}

			if (igual_pos == true) { continue; }
			else
			{
				ac.chaves_anteriores();

				//atribuir nova chave
				ac.chi[pos1ch][pos2ch] = ps.noi[pos_vizinhas[e]];
				ac.chf[pos1ch][pos2ch] = ps.nof[pos_vizinhas[e]];

				//atribuir nova posicao da chave
				ac.posicaochaves[pos1ch][pos2ch] = pos_vizinhas[e];

				//calcular FO
				soluc = ac.calculo_funcao_objetivo(pos1ch);

				//analisar caso
				if (soluc < mainsoluc)
				{
					mainsoluc = soluc;
				}
			}
		}

		pos_vizinhas.clear();
	}


	identch.clear();

	return mainsoluc;
}

float VND::v2_VND(int ch2, float incumbentv2)
{
	//troca as chaves para o vizinho do vizinho, ou seja, adjacente do adjacente
	//chavesv2 é o numero das chaves

	int identf = 0;
	vector <int> pos_vizinhas;

	bool repetepos = false;
	int pos_vetorcaso_nao_melhor = 0;

	vector <float> result_parcial; //resultado parcial da funcao objetivo
	float soluc = 0.0;
	float mainsoluc2 = 0.0;

	vector<int>auxfunc; //vetor auxiliar das funcoes
	vector <vector<int>> identch; //sempre tera duas posicoes:: 1)chi - 2)chf
	vector <vector<int>> identch_aux; //sempre tera duas posicoes:: 1)chi - 2)chf - eh o vetor auxiliar usado para encontrar as chaves

	mainsoluc2 = incumbentv2;

	// 1) localizar chaves sorteadas

	identf = 0;

	for (int j = 1; j < num_AL; j++)
	{
		for (int k = 1; k < linha_dados; k++)
		{
			if (ac.chi[j][k] != 0 && ac.chf[j][k] != 0)
			{
				identf++;

				if (identf == ch2)
				{
					auxfunc.push_back(ac.chi[j][k]);
					auxfunc.push_back(ac.chf[j][k]);

					identch.push_back(auxfunc);
					auxfunc.clear();
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
		soluc = incumbentv2;
		pos_vetorcaso_nao_melhor = ac.posicaochaves[pos1ch][pos2ch];
		bool igual_pos = false;


		for (int e = 0; e < pos_vizinhas.size(); e++)
		{
			igual_pos = false;
			for (int d = 1; d < linha_dados; d++)
			{
				if (ac.posicaochaves[pos1ch][d] != 0 && ac.posicaochaves[pos1ch][d] == pos_vizinhas[e])
				{
					igual_pos = true;
				}
			}

			if (igual_pos == true) { continue; }
			else
			{
				ac.chaves_anteriores();

				//atribuir nova chave
				ac.chi[pos1ch][pos2ch] = ps.noi[pos_vizinhas[e]];
				ac.chf[pos1ch][pos2ch] = ps.nof[pos_vizinhas[e]];

				//atribuir nova posicao da chave
				ac.posicaochaves[pos1ch][pos2ch] = pos_vizinhas[e];

				//calculo FO
				soluc = ac.calculo_funcao_objetivo(pos1ch);

				//analisar caso
				if (soluc < mainsoluc2)
				{
					mainsoluc2 = soluc;
				}
			}
		}

		pos_vizinhas.clear();
		auxfunc.clear();
		identch_aux.clear();
	}


	identch.clear();

	return mainsoluc2;

}

float VND::VND_intensificacao(int ch, float sol_incumbent)
{
	//O VND tem o objetivo de intensificar a busca do algoritmo

	float vnd_incumbent = 0.0;
	float vnd_current = 0.0;

	vnd_incumbent = sol_incumbent;
	vnd_current = sol_incumbent;

inicioVND:

	vnd_current = v1_VND(ch, vnd_incumbent);

	if (vnd_current < vnd_incumbent)
	{
		vnd_incumbent = vnd_current;
		q_vnd1++;

		//cout << "_ 1-VND: " << vnd_incumbent << endl;

		goto inicioVND;
	}

	vnd_current = v2_VND(ch, vnd_incumbent);

	if (vnd_current < vnd_incumbent)
	{
		vnd_incumbent = vnd_current;
		q_vnd2++;

		//cout << "_ 2-VND: " << vnd_incumbent << endl;

		goto inicioVND;
	}

	//fim vnd
	return vnd_incumbent;
}

float RVNS::v1_RVNS(float incumbentmainv1)
{
	//DESCRICAO: selecionar 4 chaves aleatorias do sistema e fazer o VND para estas chaves

	int ch = 0;
	float incumbentv1 = 0.0;
	float resultadov1 = 0.0;

	incumbentv1 = incumbentmainv1;


	for (int i = 1; i < 5; i++)
	{
		ch = rand() % ac.numch_SIS + 1;

		resultadov1 = vnd.VND_intensificacao(ch, incumbentv1);

		if (resultadov1 < incumbentv1)
		{
			incumbentv1 = resultadov1;
		}
	}

	return incumbentv1;
}

float RVNS::v2_RVNS(float incumbentmainv2)
{
	//DESCRICAO: sorteia dois alimentadores e faz VND em suas chaves

	float resultadov2 = 0.0;
	float incumbentv2 = 0.0;
	int al1 = 0;
	int al2 = 0;
	int cont2 = 0;
	vector<int>ch2 = {};

	al1 = al2 = 0;
	while (al1 == al2)
	{
		al1 = rand() % (num_AL - 1) + 1;
		al2 = rand() % (num_AL - 1) + 1;
	}

	//identificando chaves
	cont2 = 0;
	for (int i = 1; i < num_AL; i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			if (ac.posicaochaves[i][j] != 0)
			{
				cont2++;

				if (al1 == i || al2 == i)
				{
					ch2.push_back(cont2);
				}
			}
		}
	}

	//fazendo o VND nessas chaves dos alimentadores
	incumbentv2 = incumbentmainv2;

	for (int i = 0; i < ch2.size(); i++)
	{
		resultadov2 = vnd.VND_intensificacao(ch2[i], incumbentv2);

		if (resultadov2 < incumbentv2)
		{
			incumbentv2 = resultadov2;
		}
	}

	ch2.clear();

	return incumbentv2;
}

float RVNS::v3_RVNS(float incumbentmainv3)
{
	//DESCRICAO: seleciona dois alimentadores, seleciona uma chave para manter de um deles e faz o vnd para todas as chaves desses alimentadores
	int al1 = 0;
	int al2 = 0;
	int almantemch = 0;

	int chiv3 = 0;
	int chfv3 = 0;
	int poschv3 = 0;
	int sortchv3 = 0;
	bool cond = false;

	int chiv3b = 0;
	int chfv3b = 0;
	int poschv3b = 0;
	int sortchv3b = 0;

	vector<int>ch3 = {};
	vector<int>ch3b = {};
	int cont3 = 0;

	float incumbentv3 = 0.0;
	float resultadov3 = 0.0;
	float retorna_incumbentv3 = 0.0;

	//valor de incumbente
	retorna_incumbentv3 = incumbentmainv3;

	//sorteio alimentadores
	al1 = 0;
	al2 = 0;

	while (al1 == al2)
	{
		al1 = rand() % (num_AL - 1) + 1;
		al2 = rand() % (num_AL - 1) + 1;
	}

	//pegando chaves
	cont3 = 0;
	for (int i = 1; i < num_AL; i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			if (ac.posicaochaves[i][j] != 0)
			{
				cont3++;
				if (i == al1)
				{
					ch3.push_back(cont3);
				}
				else if (i == al2)
				{
					ch3b.push_back(cont3);
				}
			}
		}
	}

	//Paro o PRIMEIRO ALIMENTADOR - al1

	//selecionando uma chave do al1 para manter
	sortchv3 = rand() % ac.numch_AL[al1] + 1;

	chiv3 = ac.chi[al1][sortchv3];
	chfv3 = ac.chf[al1][sortchv3];
	poschv3 = ac.posicaochaves[al1][sortchv3];

	gvns.chaves_anterioresVND(al1);
	gvns.fo_anteriorVND(al1);

	//sorteando demais chaves
	gvns.sorteiochaves(ac.numch_AL[al1], fxp.camadaAL[al1], ac.posicaochaves[al1], al1);

	for (int j = 1; j < linha_dados; j++)
	{
		ac.chi[al1][j] = ps.noi[ac.posicaochaves[al1][j]];
		ac.chf[al1][j] = ps.nof[ac.posicaochaves[al1][j]];
	}


	//voltando chave
	cond = false;
	for (int j = 1; j < linha_dados; j++)
	{
		if (poschv3 == ac.posicaochaves[al1][j])
		{
			cond = true; //ja tem a chave
		}
	}

	if (cond == false)
	{
		ac.chi[al1][sortchv3] = chiv3;
		ac.chf[al1][sortchv3] = chfv3;
		ac.posicaochaves[al1][sortchv3] = poschv3;
	}

	//analisando vnd
	incumbentv3 = ac.calculo_funcao_objetivo_geral();

	for (int i = 0; i < ch3.size(); i++)
	{
		resultadov3 = vnd.VND_intensificacao(ch3[i], incumbentv3);

		if (resultadov3 < incumbentv3)
		{
			incumbentv3 = resultadov3;
		}
	}

	if (incumbentv3 > retorna_incumbentv3)
	{
		gvns.volta_chaves_anterioresVND(al1);
		gvns.volta_fo_anteriorVND(al1);
	}
	else
	{
		retorna_incumbentv3 = incumbentv3;
	}

	ch3.clear();

	//Para o SEGUNDO ALIMENTADOR - al2

	gvns.chaves_anterioresVND(al2);
	gvns.fo_anteriorVND(al2);

	//selecionando duas chaves do al2 para manter
	sortchv3 = sortchv3b = 0;
	while (sortchv3 == sortchv3b)
	{
		sortchv3 = rand() % ac.numch_AL[al2] + 1;
		sortchv3b = rand() % ac.numch_AL[al2] + 1;
	}
	
	chiv3 = ac.chi[al2][sortchv3];
	chfv3 = ac.chf[al2][sortchv3];
	poschv3 = ac.posicaochaves[al2][sortchv3];

	chiv3b = ac.chi[al2][sortchv3b];
	chfv3b = ac.chf[al2][sortchv3b];
	poschv3b = ac.posicaochaves[al2][sortchv3b];

	//sorteio do segundo alimentador
	gvns.sorteiochaves(ac.numch_AL[al2], fxp.camadaAL[al2], ac.posicaochaves[al2], al2);

	for (int j = 1; j < linha_dados; j++)
	{
		ac.chi[al2][j] = ps.noi[ac.posicaochaves[al2][j]];
		ac.chf[al2][j] = ps.nof[ac.posicaochaves[al2][j]];
	}

	//voltando as duas chaves salvas

	//chave 1
	cond = false;
	for (int j = 1; j < linha_dados; j++)
	{
		if (poschv3 == ac.posicaochaves[al2][j])
		{
			cond = true; //ja tem a chave
		}
	}

	if (cond == false)
	{
		ac.chi[al2][sortchv3] = chiv3;
		ac.chf[al2][sortchv3] = chfv3;
		ac.posicaochaves[al2][sortchv3] = poschv3;
	}

	//chave 2
	cond = false;
	for (int j = 1; j < linha_dados; j++)
	{
		if (poschv3b == ac.posicaochaves[al2][j])
		{
			cond = true; //ja tem a chave
		}
	}

	if (cond == false)
	{
		ac.chi[al2][sortchv3b] = chiv3b;
		ac.chf[al2][sortchv3b] = chfv3b;
		ac.posicaochaves[al2][sortchv3b] = poschv3b;
	}

	//analisando vnd para o segundo al
	incumbentv3 = ac.calculo_funcao_objetivo_geral();

	for (int i = 0; i < ch3b.size(); i++)
	{
		resultadov3 = vnd.VND_intensificacao(ch3b[i], incumbentv3);

		if (resultadov3 < incumbentv3)
		{
			incumbentv3 = resultadov3;
		}
	}

	if (incumbentv3 > retorna_incumbentv3)
	{
		gvns.volta_chaves_anterioresVND(al2);
		gvns.volta_fo_anteriorVND(al2);
	}
	else
	{
		retorna_incumbentv3 = incumbentv3;
	}

	ch3b.clear();

	return retorna_incumbentv3;
}

float RVNS::v4_RVNS(float incumbentmainv4)
{
	//DESCRICAO: seleciona todos os alimantadores, mantem uma chave de cada um deles e faz vnd para todas elas

	int almantemch = 0;

	int chiv4 = 0;
	int chfv4 = 0;
	int poschv4 = 0;
	int sortchv4 = 0;
	bool cond = false;

	vector<int>ch4 = {};
	int cont4 = 0;

	float incumbentv4 = 0.0;
	float resultadov4 = 0.0;
	float retorna_incumbentv4 = 0.0;

	//valor de incumbente
	retorna_incumbentv4 = incumbentmainv4;

	//analise para todos os alimentadores
	for (int e = 1; e < num_AL; e++)
	{
		//pegando chaves
		cont4 = 0;
		for (int i = 1; i < num_AL; i++)
		{
			for (int j = 1; j < linha_dados; j++)
			{
				if (ac.posicaochaves[i][j] != 0)
				{
					cont4++;
					if (i == e)
					{
						ch4.push_back(cont4);
					}
				}
			}
		}

		//selecionando uma chave do e para manter
		sortchv4 = rand() % ac.numch_AL[e] + 1;

		chiv4 = ac.chi[e][sortchv4];
		chfv4 = ac.chf[e][sortchv4];
		poschv4 = ac.posicaochaves[e][sortchv4];

		gvns.chaves_anterioresVND(e);
		gvns.fo_anteriorVND(e);

		//sorteando demais chaves
		gvns.sorteiochaves(ac.numch_AL[e], fxp.camadaAL[e], ac.posicaochaves[e], e);

		for (int j = 1; j < linha_dados; j++)
		{
			ac.chi[e][j] = ps.noi[ac.posicaochaves[e][j]];
			ac.chf[e][j] = ps.nof[ac.posicaochaves[e][j]];
		}


		//voltando chave
		cond = false;
		for (int j = 1; j < linha_dados; j++)
		{
			if (poschv4 == ac.posicaochaves[e][j])
			{
				cond = true; //ja tem a chave
			}
		}

		if (cond == false)
		{
			ac.chi[e][sortchv4] = chiv4;
			ac.chf[e][sortchv4] = chfv4;
			ac.posicaochaves[e][sortchv4] = poschv4;
		}

		//analisando vnd
		incumbentv4 = ac.calculo_funcao_objetivo_geral();

		for (int i = 0; i < ch4.size(); i++)
		{
			resultadov4 = vnd.VND_intensificacao(ch4[i], incumbentv4);

			if (resultadov4 < incumbentv4)
			{
				incumbentv4 = resultadov4;
			}
		}

		if (incumbentv4 > retorna_incumbentv4)
		{
			gvns.volta_chaves_anterioresVND(e);
			gvns.volta_fo_anteriorVND(e);
		}
		else
		{
			retorna_incumbentv4 = incumbentv4;
		}

		ch4.clear();

	}


	return retorna_incumbentv4;
}

//############################################################################################

int main()
{
	srand(static_cast <unsigned int> (time(NULL)));	//faz a aleatoriedade com base no relogio
	//srand(time(NULL));

	int itGVNS = 0;
	int simulacao = 0;

inicio_alg:

	simulacao++;

	cout << "Simulacao: " << simulacao << endl;
	cout << "\n";

	//variaveis a serem analisadas
	float current_solution = 0.0;
	float incumbent_solution = 0.0;

	rvns.q_rvns1 = 0;
	rvns.q_rvns2 = 0;
	rvns.q_rvns3 = 0;
	rvns.q_rvns4 = 0;
	vnd.q_vnd1 = 0;
	vnd.q_vnd2 = 0;

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

	//suposicao para o alimentador 1: TESTE
	/*
	ac.posicaochaves[1][1] = 7;
	ac.posicaochaves[1][2] = 9;
	ac.posicaochaves[1][3] = 14;
	*/
	

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

	//calcular valor da FO total

	incumbent_solution = ac.calculo_funcao_objetivo_geral();
	
	ac.chaves_anteriores();
	
	current_solution = incumbent_solution;

	cout << "Solucao Inicial: " << incumbent_solution << endl;
	cout << "\n";

	itGVNS = 0;

nmgvns:

	itGVNS++;
	cout << "Processo de parada - " << itGVNS << "/" << criterio_parada << endl;
	cout << "\n";

metaheuristicGVNS:

	cout << ". \n";
	current_solution = rvns.v1_RVNS(incumbent_solution);

	if (current_solution < incumbent_solution)
	{
		incumbent_solution = current_solution;
		rvns.q_rvns1++;

		cout << "1-RVNS: " << incumbent_solution << endl;

		goto metaheuristicGVNS;
	}

	cout << ". \n";
	current_solution = rvns.v2_RVNS(incumbent_solution);

	if (current_solution < incumbent_solution)
	{
		incumbent_solution = current_solution;
		rvns.q_rvns2++;

		cout << "2-RVNS: " << incumbent_solution << endl;

		goto metaheuristicGVNS;
	}

	cout << ". \n";
	current_solution = rvns.v3_RVNS(incumbent_solution);

	if (current_solution < incumbent_solution)
	{
		incumbent_solution = current_solution;
		rvns.q_rvns3++;

		cout << "3-RVNS: " << incumbent_solution << endl;

		goto metaheuristicGVNS;
	}

	cout << ". \n";
	current_solution = rvns.v4_RVNS(incumbent_solution);

	if (current_solution < incumbent_solution)
	{
		incumbent_solution = current_solution;
		rvns.q_rvns4++;

		cout << "4-RVNS: " << incumbent_solution << endl;

		goto metaheuristicGVNS;
	}

	if (itGVNS < criterio_parada)
	{
		cout << "\n";
		goto nmgvns;
	}

	////////////////////////////////////////////////////////////
	//Fim GVNS

	cout << "\n";

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
	cout << "rvns 1: " << rvns.q_rvns1 << endl;
	cout << "rvns 2: " << rvns.q_rvns2 << endl;
	cout << "rvns 3: " << rvns.q_rvns3 << endl;
	cout << "rvns 4: " << rvns.q_rvns4 << endl;
	cout << "vnd 1: " << vnd.q_vnd1 << endl;
	cout << "vnd 2: " << vnd.q_vnd2 << endl;
	cout << "\n";

	//fluxo de potencia:
	cout << "Numero de fluxo de potencia da simulacao: " << fxp.contadorFXP << endl;
	cout << "\n";

	//Solucao final
	cout << "Solucao otima: " << incumbent_solution << endl;
	cout << "\n\n";

	cout << "##########################################################################################################################" << endl;
	cout << "\n\n";

	if (simulacao < numero_simulacoes) { goto inicio_alg; }

	cout << "\n";
	cout << "FIM" << endl;

	return 0;
}