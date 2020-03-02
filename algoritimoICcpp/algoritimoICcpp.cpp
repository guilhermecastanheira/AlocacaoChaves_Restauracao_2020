//IC 2020
//ALOCAÇÃO DE CHAVES DE MANOBRA PARA RESTAURAÇÃO

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <time.h>
#include <cmath>
#include <complex>
#include <vector>
#include <string>

// DADOS DO SISTEMA DE 136 BARRAS ---------------------------------------

//Dados para o arquivo txt
#define linha_dados 157 //numero de linhas da matriz de dados +1
#define coluna_dados 9 //numero de colunas da matriz de dados +1

//Caracteristicas -----------------------------------------------------

#define num_AL 9 //numero de alimentadores +1 por conta do cpp
#define estado_restaurativo_pu 0.93 //tensao minima no estado restaurativo

//Chave a cada quantos kW?
#define parametroCH_kW 1000;

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
std::complex <float> criterio_conv = 1 * pow(10, -5);
float epsilon = abs(criterio_conv);

// Caracteristicas dos Alimentadores e Subestacoes ---------------------

int alimentadores[num_AL] = { 0, 1000, 1001, 1002, 1003, 1004, 1005, 1006, 1007 };
float capSE[num_AL] = { 0, sref, sref, sref, sref, sref, sref, sref, sref };

//NOTAS:
/*
- chaves no estado aberto eh considerado 0
- chaves fechadas é considerado 1


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

	int cadidato_aloc[linha_dados]; //candidato a alocação de chaves
	int estado_swt[linha_dados]; //estado da chave

	float dist_no[linha_dados]; //distancia entre nós

	float total_ativa = 0;
	std::complex <float> total_complexa = std::complex <float>(float (0.0), float (0.0));

	std::complex <float> lt[linha_dados]; //linha de transmissao entre nós
	std::complex <float> s_nof[linha_dados]; //potencia complexa do nof

	std::complex <float> pu_lt[linha_dados]; //linha de transmissao entre nós em pu
	std::complex <float> pu_s_nof[linha_dados]; //potencia complexa do nof em pu

	void leitura_parametros();
	void somatorio_potencia();

}ps;

class FluxoPotencia
{
public:

	int camadaAL[num_AL][linha_dados][linha_dados];
	
	int conexao_predef[linha_dados][3];

	std::complex <float> tensao_inicial = std::complex <float>(float(tensao_inicial_nos / vref), float(0)); // tensao complexa nos nós na 1 iteraçao do fluxo de potencia

	std::complex <float> corrente_pu[linha_dados];
	std::complex <float> tensao_pu[linha_dados];

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
	int posicaochaves[num_AL][linha_dados]; //vetor com as posicoes das chaves
	int adjacente_chaves[num_AL][linha_dados][linha_dados];
	int secoes_chaves[num_AL][linha_dados][linha_dados];
	int chi[num_AL][linha_dados];
	int chf[num_AL][linha_dados];


	void criterio_numero_de_chaves();
	void secoes_alimentador();
	void calculo_funcao_objetivo();
	
	

private:

	int contagem_criterio(int camada[linha_dados][linha_dados]); //criterio para a contagem de quantas chaves alocar em cada alimentador do sistema teste
	void adjacentes(int posicao[linha_dados], int adj[linha_dados][linha_dados], int alimentador); //calcula os adjacentes das chaves e da secao do alimentador
	float energia_nao_suprida(int bar_aliment[linha_dados]); //aqui se calcula a energia nao suprida para o calculo da funcao objetivo e tambem calcula a capacidade da subestacao e as condicoes de estado restaurativo
	float FO(float potencia_secao, float comprimeto_secao, float ens_isolacao);

}ac;

class GVNS
{
public:
	
	void primeiraaloc();

private:

	void sorteiochaves(int numch, int camada[linha_dados][linha_dados], int posicao_camada[linha_dados], int alimentador);
}gvns;

class RVNS:public GVNS
{
public:
	int k = 2; //numero de vizinhança do RVNS

}rvns;

class VND:public GVNS
{
public:
	int l = 2; //numero de vizinhança do VND


}vnd;


//------------------------------------------------------------

void ParametrosSistema::leitura_parametros()
{
	FILE* arquivo;

	if ((arquivo = fopen("dados136.txt", "r")) == NULL)
		return;

	for (int i = 1; i < linha_dados; i++)
	{
		fscanf(arquivo, "%d%d%f%f%f%f%d%d%f", &ps.noi[i], &ps.nof[i], &ps.lt_r[i], &ps.lt_x[i], &ps.s_nofr[i], &ps.s_nofq[i], &ps.cadidato_aloc[i], &ps.estado_swt[i], &ps.dist_no[i]);
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
	x = 0;
	for (int i = 1; i < linha_dados; i++)
	{
		if (ps.cadidato_aloc[i] == 0)
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

	std::complex <float> unit = fxp.tensao_inicial; //complexo 

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
	std::cout << "Fluxo de Potencia \n";
	std::cout << "No:" << "\t";
	std::cout << "V:" << "\t";
	std::cout << "ang:" << "\t";
	std::cout << "\t";
	std::cout << "I:" << "\t";
	std::cout << "ang:" << "\t";
	std::cout << "\n";

	for (int i = 1; i < linha_dados; i++)
	{
		std::cout << ps.nof[i] << "\t";
		std::cout << std::abs(fxp.tensao_pu[i] * vref) << "\t";
		std::cout << std::arg(fxp.tensao_pu[i]) * 180 / 3.141592 << "\t";
		//std::cout << "\t";
		std::cout << std::abs(fxp.corrente_pu[i] * iref) << "\t";
		std::cout << std::arg(fxp.corrente_pu[i]) * 180 / 3.141592 << "\n";
	}
}

void FluxoPotencia::fluxo_potencia() //alterar conforme o numero de alimentadores, modificando quantas vezes cada funcao eh chamada
{
	bool criterio_satisfeito = false;

	std::complex <float> tensao_aux[linha_dados];
	std::complex <float> corrente_aux[linha_dados];

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

		if (som == 0)
		{
			criterio_satisfeito = true;
		}
	}

} 

int AlocacaoChaves::contagem_criterio(int camada[linha_dados][linha_dados])
{
	int num_crit = 0;
	float num = 0;
	float potencia = 0;

	potencia = 0;

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
					potencia += ps.s_nofr[k];
				}
			}

		}
	}

	num = potencia / parametroCH_kW;

	num_crit = round(num);

	if (num_crit < 2) { num_crit = 2; }

	//encontando o numero estipulado
	return(num_crit);
}

void AlocacaoChaves::criterio_numero_de_chaves() //alterar conforme o numero de alimentadores, modificando quantas vezes cada funcao eh chamada
{
	for (int i = 1; i < num_AL; i++)
	{
		ac.numch_AL[i] = contagem_criterio(fxp.camadaAL[i]);
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
						if (ps.noi[k] == adj[contline][aux] && ps.cadidato_aloc[k] == 1)
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
					if (ps.noi[k] == adj[contline][aux] && ps.cadidato_aloc[k] == 1)
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

		for (int j = 2; j < linha_dados; j++)
		{
			//o j=2 porque a primeira vai ser a secao do alimentador

			for (int k = 1; k < linha_dados; k++)
			{
				//
				for (int p = 1; p < linha_dados; p++)
				{
					if (ac.secoes_chaves[i][j - 1][k] == ac.secoes_chaves[i][j][p])
					{
						ac.secoes_chaves[i][j - 1][k] = 0;
					}
				}
			}
		}
	}
}

float AlocacaoChaves::energia_nao_suprida(int bar_aliment[linha_dados])
{
	//aqui se calcula a energia nao suprida para o calculo da funcao objetivo e tambem calcula a capacidade da subestacao e as condicoes de estado restaurativo

	float potencia = 0;
	float potencia_nsup = 0;
	std::complex <float> capacidadeSE = std::complex <float>(0, 0);
	std::vector <int> posicoes_menor_estREST;
	std::string analise = "dentro do limite";

	potencia = 0.0;
	potencia_nsup = 0.0;
	analise = "dentro do limite";
	posicoes_menor_estREST.clear();

analise_potencia_nao_suprida: //aqui começa a se analisar a potencia nao suprida

	fxp.fluxo_potencia(); //separa as camadas e faz o fluxo novamente

	//analisa se existe algo fora
	for (int y = 1; y < linha_dados; y++)
	{
		if (std::abs(fxp.tensao_pu[y]) < estado_restaurativo_pu)
		{
			analise = "fora do limite";
			posicoes_menor_estREST.push_back(y);	
		}
	}

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
					if (fxp.camadaAL[i][j][k] == ps.nof[t] && fxp.camadaAL[i][j][k] != 0 && fxp.tensao_pu[t] != fxp.tensao_inicial)
					{
						potencia += ps.s_nofr[t];
						capacidadeSE += ps.s_nof[t];
					}
				}
			}
		}

		if (std::abs(capacidadeSE) > capSE[i])
		{
			analise = "fora do limite";
		}
	}

	if (analise == "dentro do limite")
	{
		potencia_nsup = ps.total_ativa - potencia;
	}
	else
	{
		//analisar qual chave pode ser aberta para diminuir a ENS
		for (int i = 0; i < posicoes_menor_estREST.size(); i++)
		{
			for (int j = 1; j < num_AL; j++)
			{
				for (int k = 1; k < linha_dados; k++)
				{
					if (ps.nof[posicoes_menor_estREST[i]] == ac.chf[j][k] || ps.nof[posicoes_menor_estREST[i]] == ac.chi[j][k] || ps.noi[posicoes_menor_estREST[i]] == ac.chf[j][k] || ps.noi[posicoes_menor_estREST[i]] == ac.chi[j][k])
					{
						if (ps.estado_swt[posicoes_menor_estREST[i]] != 0)
						{
							ps.estado_swt[posicoes_menor_estREST[i]] = 0; //fecha a chave de manobra

							goto analise_potencia_nao_suprida;//assim, temos que voltar ao inicio para nova analise	
						}
					}
				}
				
			}
		}
		
		//Ao executar este passo, não tem mais jeito... se nao estiver nas condiçoes, deve-se desligar o alimentador da falta, isso implica em somar todas as potencias
		for (int i = 1; i < linha_dados; i++)
		{
			for (int j = 1; j < linha_dados; j++)
			{
				if (bar_aliment[i] == ps.nof[j])
				{
					potencia_nsup += ps.nof[j];
				}
			}
		}
	}

	posicoes_menor_estREST.clear();

	return(potencia_nsup);

}

float AlocacaoChaves::FO(float potencia_secao, float comprimento, float ens_isolacao )
{
	float resultado = 0.0;

	resultado = (taxa_falhas * comprimento * (potencia_secao * custoKWh * tempo_falha)) + (ens_isolacao*tempo_isolacao);

	return(resultado);
}

void AlocacaoChaves::calculo_funcao_objetivo()
{
	float comprimento_secao = 0.0;
	float potencia_W = 0.0;
	float valorFO = 0.0;
	float chamadaFO = 0.0;
	float potencia_isolacao = 0.0;

	bool condicaoFOR = true;
	std::vector <int> posicao;
	std::vector <float> potencia_nao_suprida;
	
	posicao.clear();
	potencia_nao_suprida.clear();

	fxp.conexao_alimentadores(); //define as conexoes pre-existentes

	//deve-se analisar todas as secoes

	//alimentador i
	for (int i = 1; i < num_AL; i++)
	{
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
					}
				}
			}
			
			//potencia nao suprida: nesta parte é necessario analisar se o alimentador adjacente consegue suprir a energia das seçoes adjacentes a chave
			
				// 1) isolar falha
			for (int k = 1; k < linha_dados; k++)
			{
				for (int y = 1; y < linha_dados; y++)
				{
					if (ac.secoes_chaves[i][j][k] == ac.chf[i][y] || ac.secoes_chaves[i][j][k] == ac.chi[i][y])
					{
						ps.estado_swt[ac.posicaochaves[i][y]] = 0;
					}
				}
			}
			
			for (int k = 1; k < linha_dados; k++)
			{
				for (int t = 1; t < linha_dados; t++)
				{
					if (ac.secoes_chaves[i][j][k] == ps.nof[t])
					{
						ps.estado_swt[t] = 0;
					}
				}
			}

				// 2) fechar chaves da reconfiguracao para adicionar cargas ao alimentador adjacente
			for (int k = 1; k < linha_dados; k++)
			{
				for (int y = 1; y < linha_dados; y++)
				{
					if (ac.adjacente_chaves[i][j][k] == ps.nof[y] && ps.cadidato_aloc[y] == 0)
					{
						posicao.push_back(y);
					}

					if (ac.adjacente_chaves[i][j][k] == ps.noi[y] && ps.cadidato_aloc[y] == 0)
					{
						posicao.push_back(y);
					}
				}
			}

				// 3) fazer o fluxo de potencia, para todas as condições das chaves de remanejamento de cargas
			for (int y = 0; y < posicao.size(); y++)
			{
				ps.estado_swt[posicao[y]] = 1;

				//porem existem excessoes:
				for (int k = 1; k < linha_dados; k++)
				{
					if (ac.secoes_chaves[i][j][k] == ps.nof[posicao[y]] || ac.secoes_chaves[i][j][k] == ps.nof[posicao[y]])
					{
						ps.estado_swt[posicao[y]] = 0;
					}
				}

				//rodar o fluxo de potencia
				//fxp.fluxo_potencia();

				//calcular a energia nao suprida
				potencia_W = energia_nao_suprida(ac.adjacente_chaves[i][1]);

				//jogar no vetor de comparacao
				potencia_nao_suprida.push_back(potencia_W);

				//voltar a chave de remanejamento para estado aberto
				ps.estado_swt[posicao[y]] = 0;
			}

			//apos analisar todas as possiveis ligacoes, pega-se a de melhor resultado

			//4) E se não houver remanejamento??? tera que desligar toda a secao
			if (posicao.size() == 0)
			{
				for (int k = 1; k < linha_dados; k++)
				{
					for (int y = 1; y < linha_dados; y++)
					{
						if (ac.secoes_chaves[i][j][k] == ps.nof[y])
						{
							potencia_W = potencia_W + ps.s_nofr[y];
						}
					}
				}
			}
			else
			{
				//tem chave para remanejamento
				potencia_W = potencia_nao_suprida[0]; //primeiro valor encontrado

				//realizando a comparacao para ver qual o maior
				for (int y = 0; y < potencia_nao_suprida.size(); y++)
				{
					if (potencia_nao_suprida[y] <= potencia_W)
					{
						potencia_W = potencia_nao_suprida[y];
					}
				}
			}
			
			//calculando a potencia da seção para isolar
			for (int k = 1; k < linha_dados; k++)
			{
				for (int y = 1; y < linha_dados; y++)
				{
					if (ac.adjacente_chaves[i][1][k] == ps.nof[y]) //pegar a subestação inteira
					{
						potencia_isolacao = potencia_isolacao + ps.s_nofr[y];
					}
				}
			}

			//se a potencia de isolação for igual a potencia_W, significa que nao adianta fazer as manobras, sendo assim, zerar a potencia de isolacao
			if (potencia_W == potencia_isolacao)
			{
				potencia_isolacao = 0.0;
			}

			//limpar vetores comparadores
			posicao.clear();
			potencia_nao_suprida.clear();

			//deve-se calcular a funcao objetivo com os resultados obtidos
			chamadaFO = ac.FO(potencia_W, comprimento_secao, potencia_isolacao);
			valorFO = valorFO + chamadaFO;

			//se deve ler os parametros novamente para nao zerar tudo
			ps.leitura_parametros();
		}
	}

	std::cout <<"FO: " << valorFO << std::endl;
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
				if (sort == ps.nof[j] && ps.cadidato_aloc[j] == 1 && ps.noi[j] != alimentador)
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

//############################################################################################

int main()
{

	srand(static_cast <unsigned int> (time(NULL)));	//faz a aleatoriedade com base no relogio

	//faz a leitura dos parametros do circuito
	ps.leitura_parametros();

	//somatorio da potencia total do sistema
	ps.somatorio_potencia();

	//resolve o fluxo de potencia
	fxp.fluxo_potencia();

	//tirando de pu, para conferir
	fxp.valores_nominais_tensao();

	//define quantas chaves serao alocadas em cada alimentador
	ac.criterio_numero_de_chaves();

	//primeira alocacao: esta eh feita de forma aleatoria
	gvns.primeiraaloc();

	//inicio do GVNS

	//primeiras chaves
	for (int i = 1; i < num_AL; i++)
	{
		for (int j = 1; j < linha_dados; j++)
		{
			ac.chi[i][j] = ps.noi[ac.posicaochaves[i][j]];
			ac.chf[i][j] = ps.nof[ac.posicaochaves[i][j]];
		}
	}

	ac.secoes_alimentador();

	ac.calculo_funcao_objetivo(); //last edit


	// last update: 01/03/2020 - começar a fazer as vizinhanças

}