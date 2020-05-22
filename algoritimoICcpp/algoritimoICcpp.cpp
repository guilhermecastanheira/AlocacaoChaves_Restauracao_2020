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
float tensao_inicial_nos = vref * 1.05;

//critério de convergencia do fluxo de potencia
complex <float> criterio_conv = 1 * pow(10, -7);
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

	int noi[linha_dados]; //nó inicial
	int nof[linha_dados]; //nó final

	float lt_r[linha_dados]; //parte real da lt
	float lt_x[linha_dados]; //parte imaginaria da lt

	float s_nofr[linha_dados]; //potencia complexa do nof real
	float s_nofq[linha_dados]; //potencia complexa do nof img

	int candidato_aloc[linha_dados]; //candidato a alocação de chaves
	int estado_swt[linha_dados]; //estado da chave

	float dist_no[linha_dados]; //distancia entre nós

	float potencia_al[num_AL]; //potencia de cada alimentador

	float total_ativa = 0;
	complex <float> total_complexa = complex <float>(float (0.0), float (0.0));

	complex <float> lt[linha_dados]; //linha de transmissao entre nós
	complex <float> s_nof[linha_dados]; //potencia complexa do nof

	complex <float> pu_lt[linha_dados]; //linha de transmissao entre nós em pu
	complex <float> pu_s_nof[linha_dados]; //potencia complexa do nof em pu

	void leitura_parametros();
	void somatorio_potencia();

}ps;

class FluxoPotencia
{
public:

	int contadorFXP = 0; //conta quantas vzs realizou o processo de fluxo de potencia

	int camadaAL[num_AL][linha_dados][linha_dados];
	
	int conexao_predef[linha_dados][3];

	complex <float> tensao_inicial = complex <float>(float(tensao_inicial_nos / vref), float(0)); // tensao complexa nos nós na 1 iteraçao do fluxo de potencia

	complex <float> corrente_pu[linha_dados];
	complex <float> tensao_pu[linha_dados];

	void valores_nominais_tensao();
	void conexao_alimentadores();
	void fluxo_potencia();


private:

	void camadas(int alimentador, int camadaalimentador[linha_dados][linha_dados]);
	void backward_sweep(int camadaAL[linha_dados][linha_dados]);
	void forward_sweep(int alimentador, int camada[linha_dados][linha_dados]);

}fxp;

class AlocacaoChaves
{
public:

	int numch_AL[num_AL]; //numero de chaves por alimentador seguindo o criterio estipulado
	int numch_SIS = 0;
	int posicaochaves[num_AL][linha_dados]; //vetor com as posicoes das chaves
	int adjacente_chaves[num_AL][linha_dados][linha_dados];
	int secoes_chaves[num_AL][linha_dados][linha_dados];
	int chi[num_AL][linha_dados]; //barra inicial da chave
	int chf[num_AL][linha_dados]; //barra final da chave
	int antchi[num_AL][linha_dados]; //barra inicial da chave
	int antchf[num_AL][linha_dados]; //barra final da chave

	void chaves_anteriores();
	void volta_chaves_anteriores();
	void criterio_numero_de_chaves();
	void secoes_alimentador();
	float calculo_funcao_objetivo();

	
	

private:

	tuple <int, float> contagem_criterio(int camada[linha_dados][linha_dados]); //criterio para a contagem de quantas chaves alocar em cada alimentador do sistema teste
	void adjacentes(int posicao[linha_dados], int adj[linha_dados][linha_dados], int alimentador); //calcula os adjacentes das chaves e da secao do alimentador
	float energia_nao_suprida(int bar_aliment[linha_dados]); //aqui se calcula a energia nao suprida para o calculo da funcao objetivo e tambem calcula a capacidade da subestacao e as condicoes de estado restaurativo
	float FO(float potencia_secao, float comprimeto_secao, float ens_isolacao);
	float energia_suprida(float potencia_al_original, int AL);

}ac;

class GVNS
{
public:

	float current_solution = 0.0;
	float incumbent_solution = 0.0;

	float vnd_current = 0.0;
	float vnd_incumbent = 0.0;
	
	void primeiraaloc();
	float v1_RVNS(); //sortear duas chaves e aplicar VND
	float v2_RVNS(); //escolher todas as chaves de um alimentador sorteado e aplicar VND
	float v3_RVNS(); //sortear duas chaves quaisquer para outra posição e aplicar VND
	float v4_RVNS(); //sortear 2 < n < numero_max_chaves_sistema e reposiciona-las para outra posição, e aplicar VND nas n chaves 

	float v1_VND(vector<int>chavesv1); //mover para adjacente
	float v2_VND(vector<int>chavesv2); //mover para adjacente do adjacente

	float VND(vector<int>chaves);

private:

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

void FluxoPotencia::conexao_alimentadores()
{
	// ramos pre existentes no sistema, seria as linhas tracejadas no sistema
	int x = 0;
	x = 1;
	for (int i = 1; i < linha_dados; i++)
	{
		if (ps.candidato_aloc[i] == 0)
		{
			fxp.conexao_predef[x][1] = ps.noi[i];
			fxp.conexao_predef[x][2] = ps.nof[i];
			x++;
		}
	}
}

void FluxoPotencia::backward_sweep(int camadaAL[linha_dados][linha_dados])
{
	//atribuição das correntes

	for (int i = 1; i < linha_dados; i++)
	{
		for (int k = 1; k < linha_dados; k++)
		{
			for (int j = 1; j < linha_dados; j++)
			{
				if (ps.nof[i] == camadaAL[k][j] && camadaAL[k][j] != 0 && ps.estado_swt[i] == 1)
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
			for (int j = 1; j < linha_dados; j++)
			{
				if (camadaAL[i][k] == ps.nof[j] && camadaAL[i][k] != 0 && ps.estado_swt[j] == 1)
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
		cout << arg(fxp.tensao_pu[i]) * 180 / 3.141592 << "\n";
		
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

	//definir as camadas de fxp.camadas de cada alimentador

	for (int i = 1; i < num_AL; i++)
	{
		camadas(alimentadores[i], fxp.camadaAL[i]);
	}

	//tensoes iniciais nas barras

	for (int i = 1; i < linha_dados; i++)
	{
		tensao_pu[i] = fxp.tensao_inicial;
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
			convergencia_tensao[k] = abs(fxp.tensao_pu[k] - tensao_aux[k]);
			convergencia_corrente[k] = abs(fxp.corrente_pu[k] - corrente_aux[k]);
		}

		for (int k = 1; k < linha_dados; k++)
		{
			if (convergencia_tensao[k] >= epsilon)
			{
				som += 1;
			}

			if (convergencia_corrente[k] >= epsilon)
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

float AlocacaoChaves::energia_suprida(float potencia_al_original, int AL)
{
	bool analise;
	float energia_sup = 0.0;
	complex <float> capacidadeSE;
	float energia_n_sup = 0.0;
	float sum_pot = 0.0;
	
	fxp.fluxo_potencia();

	// ANALISANDO FLUXO DE POTENCIA:

	analise = true; //verdadeiro se atende as condiçoes, falso caso contrario

	// 1. Níveis de tensao devem estar acima de 0.93pu
	for (int y = 1; y < linha_dados; y++)
	{
		if (abs(fxp.tensao_pu[y]) < estado_restaurativo_pu)
		{
			analise = false;
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
						sum_pot += ps.s_nofr[t]; //aqui é apenas um somador de potencias
						capacidadeSE += ps.s_nof[t];
					}
				}
			}
		}

		if (abs(capacidadeSE) > capSE[i])
		{
			analise = false;
		}
	}

	if (analise == true)
	{
		energia_n_sup = ps.total_ativa - sum_pot;
		energia_sup = potencia_al_original - energia_n_sup; //energia suprida pela chave de remanejamento
		return(energia_sup);
	}
	else
	{
		energia_sup = 0.0;
		return(energia_sup);
	}
}

float AlocacaoChaves::FO(float potencia_secao, float comprimento, float ens_isolacao )
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

	vector<int>::iterator itr_s1;
	vector<int>::iterator itr_s2;

	bool condicaoFOR = true;

	vector <int> secao;
	vector <int> posicao;
	vector <vector<int>> analise_remanejamento;
	vector <vector<int>> remanej_cargas;
	
	posicao.clear();
	secao.clear();
	analise_remanejamento.clear();
	remanej_cargas.clear();

	//deve-se analisar todas as secoes para os valores da funcao obj

	//alimentador i
	for (int i = 1; i < num_AL; i++)
	{
		ps.leitura_parametros();

		//secao j
		for (int j = 1; j < linha_dados; j++)
		{
			//ver se vale a pena fazer o laço, se o vetor estiver zerado é só custo computacional a toa
			condicaoFOR = false;

			for (int k = 0; k < linha_dados; k++)
			{
				if (ac.secoes_chaves[i][j][k] != 0)
				{
					condicaoFOR = true;
					break;
				}
			}

			if (condicaoFOR == false) { continue; }

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

				if(cont_nao0 == 0)
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

						potencia_W = ac.energia_suprida(ps.potencia_al[i],i);

						ps.estado_swt[coluna] = 0;

						if (potencia_W != 0)
						{
							break;
						}
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

			chamadaFO = FO(ens, comprimento_secao, potencia_isolacao);

			valorFO += chamadaFO;
	
			ps.leitura_parametros();
			fxp.fluxo_potencia();
		}
	}

	cout <<"FO: " << valorFO << endl;

	return(valorFO);
}

void AlocacaoChaves::chaves_anteriores()
{
	for (int i = 1; i < num_AL; i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			antchi[i][j] = chi[i][j];
			antchf[i][j] = chf[i][j];
		}
	}
}

void AlocacaoChaves::volta_chaves_anteriores()
{
	for (int i = 1; i < num_AL; i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			chi[i][j] = antchi[i][j];
			chf[i][j] = antchf[i][j];
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

	vector <float> result_parcial; //resultado parcial da funcao objetivo
	float soluc = 0.0;

	vector<int>auxfunc; //vetor auxiliar das funcoes
	vector <vector<int>> identch; //sempre tera duas posicoes:: 1)chi - 2)chf


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

		ac.chaves_anteriores(); //salva chaves anteriores

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
		int pos_vetorcaso_nao_melhor = ac.posicaochaves[pos1ch][pos2ch];

		for (int t = 0; t < pos_vizinhas.size(); t++)
		{
			//atribuir nova chave
			ac.chi[pos1ch][pos2ch] = ps.noi[pos_vizinhas[t]];
			ac.chf[pos1ch][pos2ch] = ps.nof[pos_vizinhas[t]];

			//atribuir nova posicao da chave
			ac.posicaochaves[pos1ch][pos2ch] = pos_vizinhas[t];

			//refazer novas camadas
			ac.secoes_alimentador();

			soluc = ac.calculo_funcao_objetivo();
			result_parcial.push_back(soluc);
		}

		//pegar menor solucao

		soluc = result_parcial[0];
		int pos = 0;

		for (int t = 0; t < result_parcial.size(); t++)
		{
			if (soluc > result_parcial[t])
			{
				soluc = result_parcial[t];
				pos = t;
			}
		}

		// 4) comparar e trocar chave se necessario

		if (soluc < gvns.vnd_current)
		{
			gvns.vnd_current = soluc;

			ac.chi[pos1ch][pos2ch] = ps.noi[pos_vizinhas[pos]];
			ac.chf[pos1ch][pos2ch] = ps.nof[pos_vizinhas[pos]];

			ac.posicaochaves[pos1ch][pos2ch] = pos_vizinhas[pos];

		}
		else
		{
			ac.posicaochaves[pos1ch][pos2ch] = pos_vetorcaso_nao_melhor;
			ac.volta_chaves_anteriores();

		}


		//por fim, limpar vetores
		result_parcial.clear();
		pos_vizinhas.clear();
	}

	

	return(gvns.vnd_current);

}

float GVNS::v2_VND(vector <int> chavesv2)
{
	//troca as chaves para o vizinho do vizinho, ou seja, adjacente do adjacente
	//chavesv2 é o numero das chaves

	int identf = 0;
	vector <int> pos_vizinhas;

	bool repetepos = false;

	vector <float> result_parcial; //resultado parcial da funcao objetivo
	float soluc = 0.0;

	vector<int>auxfunc; //vetor auxiliar das funcoes
	vector <vector<int>> identch; //sempre tera duas posicoes:: 1)chi - 2)chf
	vector <vector<int>> identch_aux; //sempre tera duas posicoes:: 1)chi - 2)chf - eh o vetor auxiliar usado para encontrar as chaves


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
		///////////////////////////////////////////////////////////////////////////////

		//consertar esta parte - faz o adjacente

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
						for (int j = 0; j < identch.size(); j++)
						{
							for (int z = 0; z < identch[j].size(); z++)
							{
								if (identch[j][z] == ps.noi[k] || identch[j][z] == ps.nof[k])
								{
									partedachave = true;
								}
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
						for (int j = 0; j < identch.size(); j++)
						{
							for (int z = 0; z < identch[j].size(); z++)
							{
								if (identch[j][z] == ps.noi[k] || identch[j][z] == ps.nof[k])
								{
									partedachave = true;
								}
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
						for (int j = 0; j < identch.size(); j++)
						{
							for (int z = 0; z < identch[j].size(); z++)
							{
								if (identch[j][z] == ps.noi[k] || identch[j][z] == ps.nof[k])
								{
									partedachave = true;
								}
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
						for (int j = 0; j < identch.size(); j++)
						{
							for (int z = 0; z < identch[j].size(); z++)
							{
								if (identch[j][z] == ps.noi[k] || identch[j][z] == ps.nof[k])
								{
									partedachave = true;
								}
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

		//////////////////////////////////////////////////////////////////////////////////////

		//3) Analisar qual o melhor vizinho

		ac.chaves_anteriores(); //salva chaves anteriores

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
		int pos_vetorcaso_nao_melhor = ac.posicaochaves[pos1ch][pos2ch];

		for (int t = 0; t < pos_vizinhas.size(); t++)
		{
			//atribuir nova chave
			ac.chi[pos1ch][pos2ch] = ps.noi[pos_vizinhas[t]];
			ac.chf[pos1ch][pos2ch] = ps.nof[pos_vizinhas[t]];

			//atribuir nova posicao da chave
			ac.posicaochaves[pos1ch][pos2ch] = pos_vizinhas[t];

			//refazer novas camadas
			ac.secoes_alimentador();

			soluc = ac.calculo_funcao_objetivo();
			result_parcial.push_back(soluc);
		}

		//pegar menor solucao - consequentemente a melhor

		soluc = result_parcial[0];
		int pos = 0;

		for (int t = 0; t < result_parcial.size(); t++)
		{
			if (soluc > result_parcial[t])
			{
				soluc = result_parcial[t];
				pos = t;
			}
		}

		// 4) comparar e trocar chave se necessario

		if (soluc < gvns.vnd_current)
		{
			gvns.vnd_current = soluc;

			ac.chi[pos1ch][pos2ch] = ps.noi[pos_vizinhas[pos]];
			ac.chf[pos1ch][pos2ch] = ps.nof[pos_vizinhas[pos]];

			ac.posicaochaves[pos1ch][pos2ch] = pos_vizinhas[pos];

		}
		else
		{
			ac.posicaochaves[pos1ch][pos2ch] = pos_vetorcaso_nao_melhor;
			ac.volta_chaves_anteriores();
		}


		//por fim, limpar vetores
		result_parcial.clear();
		pos_vizinhas.clear();
		auxfunc.clear();
		identch_aux.clear();
	}



	return(gvns.vnd_current);


}

float GVNS::VND(vector <int> chaves)
{
	//O VND tem o objetivo de intensificar a busca do algoritmo

	gvns.vnd_incumbent = ac.calculo_funcao_objetivo();
	gvns.vnd_current = gvns.vnd_incumbent;

inicioVND:

	gvns.vnd_current = gvns.v1_VND(chaves);

	if (gvns.vnd_current < gvns.vnd_incumbent)
	{
		gvns.vnd_incumbent = gvns.vnd_current;
		goto inicioVND;
	}

	gvns.vnd_current = gvns.v2_VND(chaves);

	if (gvns.vnd_current < gvns.vnd_incumbent)
	{
		gvns.vnd_incumbent = gvns.vnd_current;
		goto inicioVND;
	}

	//acaba o vnd

	return(gvns.vnd_incumbent);
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
	while (sort1==sort2)
	{
		aleatorio = rand() % ac.numch_SIS + 1;
		sort1 = aleatorio; //primeiro sorteio

		aleatorio = rand() % ac.numch_SIS + 1;
		sort2 = aleatorio; //segundo sorteio
	}

	mudarch.push_back(sort1);
	mudarch.push_back(sort2);

	solution = gvns.VND(mudarch);

	mudarch.clear();

	return(solution);
	
}


//############################################################################################

int main()
{
	srand(static_cast <unsigned int> (time(NULL)));	//faz a aleatoriedade com base no relogio

	int run_alg = 0;
	int itGVNS = 0;

	vector <float> atualizacaoFO;
	vector <vector<float>> matriz_atualizacaoFO;
	vector <int> iteracoesGVNS;

rodar_dnv:

	

	//--------------------------------------------------------------------------

	//variaveis a serem analisadas
	gvns.current_solution = 0.0;
	gvns.incumbent_solution = 0.0;

	run_alg++;

	//inicializações
	ac.numch_SIS = 0;

	for (int i = 1; i < num_AL; i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			ac.chi[i][j] = 0;
			ac.chf[i][j] = 0;
		}
	}

	//faz a leitura dos parametros do circuito
	ps.leitura_parametros();

	//somatorio da potencia total do sistema
	ps.somatorio_potencia();

	//conexoes para remanejamento de cargas
	fxp.conexao_alimentadores();

	//resolve o fluxo de potencia
	fxp.fluxo_potencia();

	//tirando de pu, para conferir
	fxp.valores_nominais_tensao();

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

	ps.leitura_parametros();
	fxp.fluxo_potencia();
	ac.secoes_alimentador(); 

	//ac.calculo_funcao_objetivo(); 

metaheuristicGVNS:

	itGVNS++;

	gvns.current_solution = gvns.v1_RVNS();

	if (gvns.current_solution < gvns.incumbent_solution)
	{
		gvns.incumbent_solution = gvns.current_solution;
		atualizacaoFO.push_back(gvns.incumbent_solution);
		goto metaheuristicGVNS;
	}

	

	//end metaheuristicGVNS
	matriz_atualizacaoFO.push_back(atualizacaoFO);
	atualizacaoFO.clear();

	iteracoesGVNS.push_back(itGVNS);
	itGVNS = 0;

	if (run_alg < num_run_alg) { goto rodar_dnv; }

	//getting values


	//AQUI FAZER A ESTATISTICA DOS DADOS --- fazer depois de acabar de programar a metaheuristica





	
	cout << "\n" << "Contador Total Fluxo de Potencia: " << fxp.contadorFXP << endl;
}