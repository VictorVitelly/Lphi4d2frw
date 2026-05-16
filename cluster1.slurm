#!/bin/bash
#SBATCH --job-name=Lphi4frw                    # Nombre del trabajo
#SBATCH --output=logs/Lphi4frw.log                   # Archivo de registro de salida
#SBATCH --error=logs/Lphi4frw.err                    # Archivo de registro de errores
#SBATCH --partition=QuantPhysMC         # Nombre de la partición o cola de trabajos
#SBATCH --nodes=1                               # Número de nodos a utilizar (puedes cambiarlo)
#SBATCH --ntasks-per-node=1             # Número de tareas por nodo (1 para ejecución serial)
#SBATCH --cpus-per-task=1               # Número de CPUs por tarea (puedes cambiarlo)
#SBATCH --mem=4G                        # Memoria RAM necesaria (puedes cambiarlo)
#SBATCH --time=2-00

# Comando para ejecutar tu programa

cd ~/Lphi4d2frw

SRC=src
BIN=bin
DAT=data
TARGET=main

(echo 1) |$BIN/$TARGET
