//IC 2020
//ALOCAÇÃO DE CHAVES DE MANOBRA PARA RESTAURAÇÃO

#define _CRT_SECURE_NO_WARNINGS

#include <iostream>
#include <time.h>
#include <stdio.h>
#include <stdlib.h>
#include <cmath>
#include <complex>
#include <vector>

// DADOS DO SISTEMA ---------------------------------------------------

//Dados para o arquivo txt
#define linha_dados 157 //numero de linhas da matriz de dados +1
#define coluna_dados 9 //numero de colunas da matriz de dados +1

//Caracteristicas -----------------------------------------------------

#define num_AL 9 //numero de alimentadores +1 por conta do cpp

//alimentadores das subestações
#define a1 1000 //alimentador 1
#define a2 1001 //alimentador 2
#define a3 1002 //alimentador 3
#define a4 1003 //alimentador 4
#define a5 1004 //alimentador 5
#define a6 1005 //alimentador 6
#define a7 1006 //alimentador 7
#define a8 1007 //alimentador 8

//Chave a cada quantos kW?
#define parametroCH_kW 1000;

// dados dos condutores a serem usados (catálogo Nexans): 
std::complex <float> t_raven = std::complex <float>(0.7208, 0.4186); //T-Raven 

//Caracteristicas Fluxo de Potencia ------------------------------------

//valores base, ou, de referencia
float sref = 100 * pow(10, 6); //100MVA
float vref = 13800; //13.8kV
float zref = (vref * vref) / sref;
float iref = sref / vref;

//tensao inicial nos nós do sistema
float tensao_inicial_nos = 13800;

//critério de convergencia do fluxo de potencia
std::complex <float> criterio_conv = 1 * pow(10, -8);
float epsilon = abs(criterio_conv);



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

	std::complex <float> lt[linha_dados]; //linha de transmissao entre nós
	std::complex <float> s_nof[linha_dados]; //potencia complexa do nof

	std::complex <float> pu_lt[linha_dados]; //linha de transmissao entre nós em pu
	std::complex <float> pu_s_nof[linha_dados]; //potencia complexa do nof em pu

	void leitura_parametros();

}ps;

class FluxoPotencia
{
public:

	int camadaAL1[linha_dados][linha_dados];
	int camadaAL2[linha_dados][linha_dados];
	int camadaAL3[linha_dados][linha_dados];
	int camadaAL4[linha_dados][linha_dados];
	int camadaAL5[linha_dados][linha_dados];
	int camadaAL6[linha_dados][linha_dados];
	int camadaAL7[linha_dados][linha_dados];
	int camadaAL8[linha_dados][linha_dados];

	int conexao_predef[linha_dados][3];

	std::complex <float> tensao_inicial = std::complex <float>(float(tensao_inicial_nos / vref), float(0)); // tensao complexa nos nós na 1 iteraçao do fluxo de potencia

	std::complex <float> corrente_pu[linha_dados];
	std::complex <float> tensao_pu[linha_dados];

	void valores_nominais_tensao();
	void fluxo_potencia();


private:

	void camadas(int alimentador, int camadaalimentador[linha_dados][linha_dados]);
	void conexao_alimentadores();
	void backward_sweep(int camadaAL[linha_dados][linha_dados]);
	void forward_sweep(int alimentador, int camada[linha_dados][linha_dados]);

}fxp;

class AlocacaoChaves
{
public:

	int numch_AL[num_AL]; //numero de chaves por alimentador seguindo o criterio estipulado
	int posicaochaves[num_AL][linha_dados]; //vetor com as posicoes das chaves

private:

	void criterio_numero_de_chaves();
	void contagem_criterio(int camada[linha_dados][linha_dados], int numeroestipulado);

}ac;

class MetaheuristicaGVNS
{
public:

	//void primeiraalocacaosistema(int numch);

}gvns;

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

void FluxoPotencia::camadas(int alimentador, int camadaalimentador[linha_dados][linha_dados])
{
	//zerar camada 

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

	for (int k = 1; k < linha_dados; k++)
	{
		for (int m = 1; m < linha_dados; m++)
		{
			for (int i = 1; i < linha_dados; i++)
			{

				if (ps.noi[i] == camadaalimentador[k][m])
				{
					if (ps.estado_swt[i] == 1)
					{
						camadaalimentador[x][y] = ps.nof[i];
						y += 1;
					}

				}


			}
		}

		x += 1;
		y = 1;
	}

}

void FluxoPotencia::conexao_alimentadores()
{
	// ramos pre existentes no sistema, seria as linhas tracejadas no sistema

	for (int i = 1; i < linha_dados; i++)
	{
		if (ps.cadidato_aloc[i] == 0)
		{
			fxp.conexao_predef[i][1] = ps.noi[i];
			fxp.conexao_predef[i][2] = ps.nof[i];
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

	std::complex <float> unit = 1.0; //complexo 1|0°_

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

	camadas(a1, fxp.camadaAL1);
	camadas(a2, fxp.camadaAL2);
	camadas(a3, fxp.camadaAL3);
	camadas(a4, fxp.camadaAL4);
	camadas(a5, fxp.camadaAL5);
	camadas(a6, fxp.camadaAL6);
	camadas(a7, fxp.camadaAL7);
	camadas(a8, fxp.camadaAL8);


	conexao_alimentadores(); //define as conexoes pre-existentes

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

		backward_sweep(fxp.camadaAL1);
		backward_sweep(fxp.camadaAL2);
		backward_sweep(fxp.camadaAL3);
		backward_sweep(fxp.camadaAL4);
		backward_sweep(fxp.camadaAL5);
		backward_sweep(fxp.camadaAL6);
		backward_sweep(fxp.camadaAL7);
		backward_sweep(fxp.camadaAL8);

		//2 passo: FORWARD

		forward_sweep(a1, fxp.camadaAL1);
		forward_sweep(a2, fxp.camadaAL2);
		forward_sweep(a3, fxp.camadaAL3);
		forward_sweep(a4, fxp.camadaAL4);
		forward_sweep(a5, fxp.camadaAL5);
		forward_sweep(a6, fxp.camadaAL6);
		forward_sweep(a7, fxp.camadaAL7);
		forward_sweep(a8, fxp.camadaAL8);

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

void AlocacaoChaves::contagem_criterio(int camada[linha_dados][linha_dados], int numeroestipulado)
{
	int potencia = 0;

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
					potencia =+ ps.s_nofr[k];
				}
			}

		}
	}

	//encontando o numero estipulado
	numeroestipulado = potencia / parametroCH_kW;
}

void AlocacaoChaves::criterio_numero_de_chaves() //alterar conforme o numero de alimentadores, modificando quantas vezes cada funcao eh chamada
{
	contagem_criterio(fxp.camadaAL1, ac.numch_AL[1]);
	contagem_criterio(fxp.camadaAL2, ac.numch_AL[2]);
	contagem_criterio(fxp.camadaAL3, ac.numch_AL[3]);
	contagem_criterio(fxp.camadaAL4, ac.numch_AL[4]);
	contagem_criterio(fxp.camadaAL5, ac.numch_AL[5]);
	contagem_criterio(fxp.camadaAL6, ac.numch_AL[6]);
	contagem_criterio(fxp.camadaAL7, ac.numch_AL[7]);
	contagem_criterio(fxp.camadaAL8, ac.numch_AL[8]);
}

/*
void MetaheuristicaGVNS::primeiraalocacaosistema(int numch)
{
	//na primeira alocacao, temos as chaves sorteadas para cada alimentador

	int sorteio = 0; //variavel com o valor do sorteio
	int cont = 0;

	for (int i = 1; i < numch; i++)
	{
	dnv:
		sorteio = rand() % (linha_dados - 1) + 1;

		if (i == 1 && ps.cadidato_aloc[sorteio] != 0)
		{
			ac.posicaochaves[i] = sorteio;
		}
		else if (ps.cadidato_aloc[sorteio] != 0)
		{
			cont = 0;
			for (int k = 1; k < numch; k++)
			{
				for (int j = 1; j < numch; j++)
				{
					if (ac.posicaochaves[k] == ac.posicaochaves[j] && ac.posicaochaves[j] != 0)
					{
						cont++;
					}
				}
			}

			if (cont >= i) { goto dnv; }
			else { ac.posicaochaves[i] = sorteio; }
		}
		else { goto dnv; }
	}
}
*/

//############################################################################################

int main()
{

	srand(static_cast<unsigned int> (time(NULL)));	//faz a aleatoriedade com base no relogio

	//faz a leitura dos parametros do circuito
	ps.leitura_parametros();

	//resolve o fluxo de potencia
	fxp.fluxo_potencia();

	//tirando de pu, para conferir
	fxp.valores_nominais_tensao();

	//gvns.primeiraalocacaosistema();

}