__author__ = 'Evan'

import time
import datetime
import tkinter as Tk
from edu.unco.mast.src.ga import Algorithm, CrossoverType, SelectionType, TerminationType
import PIL

class App:

    CONTROLWIDTH = 512
    CONTROLHEIGHT = 512

    DISPLAYWIDTH = 512
    DISPLAYHEIGHT = 512

    GENENUMBERS = [8, 16, 32, 64]
    POPULATIONSIZES = [10, 100, 1000, 5000]
    ELITEPOPULATIONPERCENTS = [0, 10, 25, 50, 75, 90, 100]
    GENERATIONNUMBERS = [10, 100, 1000, 5000]
    CROSSOVERCHANCES = [10, 25, 50, 75, 90, 100]
    MUTATIONSIZES = [0, 1, 2, 3]
    MUTATIONCHANCES = [10, 25, 50, 75, 90, 100]

    def __init__(self):
        self.__algorithm = None
        self.__numGenes = 8
        self.__popSize = 100
        self.__eliteSize = 50
        self.__maxGen = 100
        self.__crossoverType = "SINGLEPOINT"
        self.__crossoverChance = 100
        self.__mutationSize = 1
        self.__mutationChance = 10
        self.__selectionType = "RANDOMELITE"
        self.__terminationType = "CONVERGENCE"
        self.__optimumX = -1
        self.__optimumY = -1
        self.__convergenceConfidence = 10
        self.__permutateAttributes = False
        self.__evaluationDictionary = {}
        self.__displayImagePath = ""

        self.__root = Tk.Tk()
        self.__root.wm_title("GeneTwo Genetic Algorithm Simulation")
        self.__buildGUI(self.__root)

    def __buildGUI(self, master):
        self.__controlFrame = Tk.Frame(master, width = self.CONTROLWIDTH, height = self.CONTROLHEIGHT)
        self.__controlFrame.grid(row = 0, column = 0)
        self.__displayFrame = Tk.Frame(master, width = self.DISPLAYWIDTH, height = self.DISPLAYHEIGHT)
        self.__displayFrame.grid(row = 0, column = 1)

        self.__configButton = Tk.Button(
            self.__controlFrame, text = "Load Config File", command = self.__buildLoadDialog, width = 20
        )
        self.__configButton.grid(row = 0, column = 0)
        self.__attributeButton = Tk.Button(
            self.__controlFrame, text = "Customize Attributes", command = self.__buildAttributeDialog, width = 20
        )
        self.__attributeButton.grid(row = 1, column = 0)
        self.__controlSpacer = Tk.Frame(self.__controlFrame, height = 250)
        self.__controlSpacer.grid(row = 2, column = 0)
        self.__executeButton = Tk.Button(
            self.__controlFrame, text = "Execute", command = self.__runAlgorithm, width = 20
        )
        self.__executeButton.grid(row = 3, column = 0)
        self.__bruteForceButton = Tk.Button(
            self.__controlFrame, text = "Brute Force", command = self.__runBruteForce, width = 20
        )
        self.__bruteForceButton.grid(row = 4, column = 0)

    def __buildNotificationDialog(self, message):
        self.__notificationDialog = Tk.Toplevel()
        self.__notificationDialog.title("Message")
        self.__notificationDialog.resizable(False, False)
        self.__notificationDialog.minsize(width = 150, height = 50)

        self.__notificationLabel = Tk.Label(self.__notificationDialog, text=message)
        self.__notificationLabel.pack()

        self.__notificationDialog.grab_set()
        self.__notificationDialog.focus()

    def __buildLoadDialog(self):
        self.__loadDialog = Tk.Toplevel()
        self.__loadDialog.title("Load GeneTwo Config File")
        self.__loadDialog.resizable(False, False)

        self.__loadFrame = Tk.Frame(self.__loadDialog)
        self.__loadFrame.grid(row = 0, column = 0)

        self.__loadField = Tk.Entry(self.__loadFrame, width = 50)
        self.__loadField.insert(0, "Enter file path...")
        self.__loadField.grid(row = 0, column = 0)

        self.__loadButton = Tk.Button(
            self.__loadFrame, text = "Load", command = self.__loadConfig
        )
        self.__loadButton.grid(row = 2, column = 0)
        self.__loadDialog.grab_set()
        self.__loadDialog.focus()

    def __loadConfig(self):
        if not self.__loadDialog is None:
            file = self.__loadField.get()
            file.replace("\\", "/")
            print(file)
            content = open(file).readlines()
            self.__parseConfig(content)
            self.__loadDialog.destroy()
            self.__buildNotificationDialog("Configuration settings successfully applied")

    def __parseConfig(self, config):
        i = 0
        while (i < len(config)):
            config[i] = config[i].replace(" ", "")
            if ("ATTRIBUTES<" in config[i]):
                while (i < len(config)):
                    i += 1
                    config[i] = config[i].replace(" ", "")
                    if (">" in config[i]):
                        break
                    words = config[i].split('=')
                    if (len(words) > 1):
                        if (words[0] == "PopulationSize"):
                            self.__popSize = int(words[1])
                        elif (words[0] == "EliteSize"):
                            self.__eliteSize = int(words[1])
                        elif (words[0] == "MaxGenerations"):
                            self.__maxGen = int(words[1])
                        elif (words[0] == "CrossoverType"):
                            self.__crossoverType = words[1]
                        elif (words[0] == "CrossoverChance"):
                            self.__crossoverChance = int(words[1])
                        elif (words[0] == "MutationSize"):
                            self.__mutationSize = int(words[1])
                        elif (words[0] == "MutationChance"):
                            self.__mutationChance = int(words[1])
                        elif (words[0] == "SelectionType"):
                            self.__selectionType = words[1]
                        elif (words[0] == "TerminationType"):
                            self.__terminationType = words[1]
                        elif (words[0] == "OptimumX"):
                            self.__optimumX = int(words[1])
                        elif (words[0] == "OptimumY"):
                            self.__optimumY = int(words[1])
                        elif (words[0] == "ConvergenceConfidence"):
                            self.__convergenceConfidence = int(words[1])
                        elif (words[0] == "NumGenes"):
                            self.__numGenes = int(words[1])
            elif ("EVALUATION<" in config[i]):
                self.__evaluationDictionary = {}
                while(i < len(config)):
                    i += 1
                    config[i] = config[i].replace(" ", "")
                    if (">" in config[i]):
                        break
                    words = config[i].split('=')
                    self.__evaluationDictionary.update({words[0]:float(words[1])})
            elif ("DISPLAY<" in config[i]):
                while(i < len(config)):
                    i += 1
                    config[i] = config[i].replace(" ", "")
                    if (">" in config[i]):
                        break
                    if (words[0] == "DisplayImage"):
                            self.__displayImagePath = words[1]
            i += 1

    def __buildAttributeDialog(self):
        self.__attributeDialog = Tk.Toplevel()
        self.__attributeDialog.title("Configure Genetic Algorithm")
        self.__attributeDialog.resizable(False, False)
        self.__attributeDialog.grid_columnconfigure(0, minsize = 250)
        self.__attributeDialog.grid_columnconfigure(1, minsize = 250)

        Tk.Label(self.__attributeDialog, text = "Number of Genes").grid(row = 0, column = 0)
        self.__numGenesField = Tk.Entry(self.__attributeDialog, width = 25)
        self.__numGenesField.insert(0, self.__numGenes)
        self.__numGenesField.grid(row = 0, column = 1)

        Tk.Label(self.__attributeDialog, text = "Total Population Size").grid(row = 1, column = 0)
        self.__popSizeField = Tk.Entry(self.__attributeDialog, width = 25)
        self.__popSizeField.insert(0, self.__popSize)
        self.__popSizeField.grid(row = 1, column = 1)

        Tk.Label(self.__attributeDialog, text = "Elite Population Size").grid(row = 2, column = 0)
        self.__eliteSizeField = Tk.Entry(self.__attributeDialog, width = 25)
        self.__eliteSizeField.insert(0, self.__eliteSize)
        self.__eliteSizeField.grid(row = 2, column = 1)

        Tk.Label(self.__attributeDialog, text = "Max Generations").grid(row = 3, column = 0)
        self.__maxGenField = Tk.Entry(self.__attributeDialog, width = 25)
        self.__maxGenField.insert(0, self.__maxGen)
        self.__maxGenField.grid(row = 3, column = 1)

        Tk.Label(self.__attributeDialog, text = "Crossover Type").grid(row = 4, column = 0)
        self.__crossoverVar = Tk.StringVar(self.__attributeDialog)
        self.__crossoverVar.set(self.__crossoverType)
        Tk.OptionMenu(self.__attributeDialog, self.__crossoverVar, "SINGLEPOINT", "DOUBLEPOINT", "DOUBLEPOINTBIASED",
                      "UNIFORM", "NONE").grid(row = 4, column = 1)

        Tk.Label(self.__attributeDialog, text = "Crossover Chance").grid(row = 5, column = 0)
        self.__crossoverField = Tk.Entry(self.__attributeDialog, width = 25)
        self.__crossoverField.insert(0, self.__crossoverChance)
        self.__crossoverField.grid(row = 5, column = 1)

        Tk.Label(self.__attributeDialog, text = "Mutation Size").grid(row = 6, column = 0)
        self.__mutationSizeField = Tk.Entry(self.__attributeDialog, width = 25)
        self.__mutationSizeField.insert(0, self.__mutationSize)
        self.__mutationSizeField.grid(row = 6, column = 1)

        Tk.Label(self.__attributeDialog, text = "Mutation Chance").grid(row = 7, column = 0)
        self.__mutationChanceField = Tk.Entry(self.__attributeDialog, width = 25)
        self.__mutationChanceField.insert(0, self.__mutationChance)
        self.__mutationChanceField.grid(row = 7, column = 1)

        Tk.Label(self.__attributeDialog, text = "Selection Type").grid(row = 8, column = 0)
        self.__selectionVar = Tk.StringVar(self.__attributeDialog)
        self.__selectionVar.set(self.__selectionType)
        Tk.OptionMenu(self.__attributeDialog, self.__selectionVar, "RANDOMELITE", "CASCADINGELITE", "ROULETTEELITE",
                      "ROULETTETOTAL").grid(row = 8, column = 1)

        Tk.Label(self.__attributeDialog, text = "Termination Type").grid(row = 9, column = 0)
        self.__terminationVar = Tk.StringVar(self.__attributeDialog)
        self.__terminationVar.set(self.__terminationType)
        Tk.OptionMenu(self.__attributeDialog, self.__terminationVar, "CONVERGENCE", "OPTIMUMCONVERGENCE",
                      "NATURAL",).grid(row = 9, column = 1)

        Tk.Label(self.__attributeDialog, text = "Optimum X").grid(row = 10, column = 0)
        self.__optimumXField = Tk.Entry(self.__attributeDialog, width = 25)
        self.__optimumXField.insert(0, self.__optimumX)
        self.__optimumXField.grid(row = 10, column = 1)

        Tk.Label(self.__attributeDialog, text = "Optimum Y").grid(row = 11, column = 0)
        self.__optimumYField = Tk.Entry(self.__attributeDialog, width = 25)
        self.__optimumYField.insert(0, self.__optimumY)
        self.__optimumYField.grid(row = 11, column = 1)

        Tk.Label(self.__attributeDialog, text = "Convergence Confidence").grid(row = 12, column = 0)
        self.__convergenceField = Tk.Entry(self.__attributeDialog, width = 25)
        self.__convergenceField.insert(0, self.__convergenceConfidence)
        self.__convergenceField.grid(row = 12, column = 1)

        Tk.Label(self.__attributeDialog, text = "Permutate Attributes").grid(row = 13, column = 0)
        self.__permutateVar = Tk.BooleanVar(self.__attributeDialog)
        self.__permutateField = Tk.Checkbutton(self.__attributeDialog, var = self.__permutateVar, onvalue = True, offvalue = False)
        if (self.__permutateAttributes):
            self.__permutateField.select()
        self.__permutateField.grid(row = 13, column = 1)

        Tk.Frame(self.__attributeDialog, height = 25).grid(row = 14, column = 0)
        Tk.Frame(self.__attributeDialog, height = 25).grid(row = 14, column = 1)

        self.__setAttributesButton = Tk.Button(
            self.__attributeDialog, text = "Confirm", command = self.__setAttributes
        )
        self.__setAttributesButton.grid(row = 15, column = 0)

        self.__cancelAttributesButton = Tk.Button(
            self.__attributeDialog, text = "Cancel", command = self.__cancelAttributes
        )
        self.__cancelAttributesButton.grid(row = 15, column = 1)

        self.__attributeDialog.grab_set()
        self.__attributeDialog.focus()

    def __setAttributes(self):
         if not self.__attributeDialog is None:
             self.__numGenes = int(self.__numGenesField.get())
             self.__popSize = int(self.__popSizeField.get())
             self.__eliteSize = int(self.__eliteSizeField.get())
             self.__maxGen = int(self.__maxGenField.get())
             self.__crossoverType = self.__crossoverVar.get()
             self.__crossoverChance = int(self.__crossoverField.get())
             self.__mutationSize = int(self.__mutationSizeField.get())
             self.__mutationChance = int(self.__mutationChanceField.get())
             self.__selectionType = self.__selectionVar.get()
             self.__terminationType = self.__terminationVar.get()
             self.__optimumX = int(self.__optimumXField.get())
             self.__optimumY = int(self.__optimumYField.get())
             self.__convergenceConfidence = int(self.__convergenceField.get())
             self.__permutateAttributes = self.__permutateVar.get()
             self.__attributeDialog.destroy()

    def __cancelAttributes(self):
         if not self.__attributeDialog is None:
             self.__attributeDialog.destroy()

    def __runAlgorithm(self):
        if (len(self.__evaluationDictionary) < 1):
            self.__buildNotificationDialog(str("The evaluation settings have not yet been defined."))
            return
        #Convert string-based attributes to their corresponding enums
        if self.__permutateAttributes:
            log = open("../logs/" + str(datetime.datetime.now().date()) + str(datetime.datetime.now().timestamp()) + ".txt", 'w')
            algoNum = 0
            for terminationEnum in [0, 2]:
                for selectionEnum in range(4):
                    for crossoverEnum in range(5):
                        for mutationSize in self.MUTATIONSIZES:
                            for crossoverChance in self.CROSSOVERCHANCES:
                                for mutationChance in self.MUTATIONCHANCES:
                                    for geneNumber in self.GENENUMBERS:
                                        for populationSize in self.POPULATIONSIZES:
                                            for elitePopulationPercent in self.ELITEPOPULATIONPERCENTS:
                                                for generationNumber in self.GENERATIONNUMBERS:
                                                    if (elitePopulationPercent == 0 and selectionEnum != 3):
                                                        continue
                                                    print("Writing to log")
                                                    log.write("Algorithm " + str(algoNum) + ": \n")
                                                    log.write("-NumGenes=" + str(geneNumber) + " \n")
                                                    log.write("-PopSize=" + str(populationSize) + " \n")
                                                    log.write("-ElitePop%" + str(elitePopulationPercent) + "\n")
                                                    log.write("-MaxGens=" + str(geneNumber) + " \n")
                                                    log.write("-CrossoverType=" + str(crossoverEnum) + " \n")
                                                    log.write("-CrossoverChance=" + str(crossoverChance) + " \n")
                                                    log.write("-MutationSize=" + str(mutationSize) + " \n")
                                                    log.write("-MutationChance=" + str(mutationChance) + " \n")
                                                    log.write("-SelectionType=" + str(selectionEnum) + "\n")
                                                    log.write("-TerminationType=" + str(terminationEnum) + "\n")
                                                    print("written")
                                                    for i in range(3):
                                                        self.__algorithm = Algorithm(geneNumber, populationSize, generationNumber,
                                                                                     crossoverEnum, crossoverChance, mutationSize,
                                                                                     mutationChance, round(populationSize * (elitePopulationPercent / 100)),
                                                                                     selectionEnum, terminationEnum, -1, -1, 10, self.__evaluationDictionary)
                                                        startTime = time.clock() * 1000
                                                        done = False
                                                        while not done:
                                                            done = self.__algorithm.nextGen()
                                                        log.write("Run " + str(i + 1) + ": \n")
                                                        log.write("Milliseconds elapsed: " + str((time.clock() * 1000) - startTime) + "\n")
                                                        log.write("Final Population: Generation " + str(self.__algorithm.getGenerationsPassed()))
                                                        evalData = self.__algorithm.getEvalData()
                                                        log.write("Average X: " + str(evalData[0]))
                                                        log.write("Average Y: " + str(evalData[1]))
                                                        log.write("Average Fitness: " + str(evalData[2]))
                                                        log.write("Best Fitness Scores: " + str(evalData[3]))
                                                    print("Algorithm",algoNum,"completed.")
                                                    algoNum += 1
        else:
            crossoverType = -1
            if (self.__crossoverType == "SINGLEPOINT"):
                crossoverType = CrossoverType.SINGLEPOINT
            elif (self.__crossoverType == "DOUBLEPOINT"):
                crossoverType = CrossoverType.DOUBLEPOINT
            elif (self.__crossoverType == "DOUBLEPOINTBIASED"):
                crossoverType = CrossoverType.DOUBLEPOINTBIASED
            elif (self.__crossoverType == "UNIFORM"):
                crossoverType = CrossoverType.UNIFORM
            elif (self.__crossoverType == "NONE"):
                crossoverType = CrossoverType.NONE
            else:
                self.__buildNotificationDialog("Crossover type " + self.__crossoverType + " is not properly defined.")
                return
            selectionType = -1
            if (self.__selectionType == "RANDOMELITE"):
                selectionType = SelectionType.RANDOMELITE
            elif (self.__selectionType == "CASCADINGELITE"):
                selectionType = SelectionType.CASCADINGELITE
            elif (self.__selectionType == "ROULETTEELITE"):
                selectionType = SelectionType.ROULETTEELITE
            elif (self.__selectionType == "ROULETTETOTAL"):
                selectionType = SelectionType.ROULETTETOTAL
            else:
                self.__buildNotificationDialog("Selection type " + self.__selectionType + " is not properly defined.")
                return
            terminationType = -1
            if (self.__terminationType == "CONVERGENCE"):
                terminationType = TerminationType.CONVERGENCE
            elif (self.__terminationType == "OPTIMUMCONVERGENCE"):
                terminationType = TerminationType.OPTIMUMCONVERGENCE
            elif (self.__terminationType == "NATURAL"):
                terminationType = TerminationType.NATURAL
            else:
                self.__buildNotificationDialog("Termination type " + self.__terminationType + " is not properly defined.")
                return
            self.__algorithm = Algorithm(self.__numGenes, self.__popSize, self.__maxGen, crossoverType, self.__crossoverChance, self.__mutationSize,
                                         self.__mutationChance, self.__eliteSize, selectionType, terminationType, self.__optimumX,
                                         self.__optimumY, self.__convergenceConfidence, self.__evaluationDictionary)
            startTime = time.clock() * 1000
            done = False
            while not done:
                done = self.__algorithm.nextGen()
            print("Milliseconds elapsed:",(time.clock() * 1000) - startTime)
            print("Final Population: Generation",self.__algorithm.getGenerationsPassed())
            print(self.__algorithm.getEvalData())

    def __runBruteForce(self):
        if (len(self.__evaluationDictionary) < 1):
            self.__buildNotificationDialog(str("The evaluation settings have not yet been defined."))
            return
        crossoverType = -1
        if (self.__crossoverType == "SINGLEPOINT"):
            crossoverType = CrossoverType.SINGLEPOINT
        elif (self.__crossoverType == "DOUBLEPOINT"):
            crossoverType = CrossoverType.DOUBLEPOINT
        elif (self.__crossoverType == "DOUBLEPOINTBIASED"):
            crossoverType = CrossoverType.DOUBLEPOINTBIASED
        elif (self.__crossoverType == "UNIFORM"):
            crossoverType = CrossoverType.UNIFORM
        elif (self.__crossoverType == "NONE"):
            crossoverType = CrossoverType.NONE
        else:
            self.__buildNotificationDialog("Crossover type " + self.__crossoverType + " is not properly defined.")
            return
        selectionType = -1
        if (self.__selectionType == "RANDOMELITE"):
            selectionType = SelectionType.RANDOMELITE
        elif (self.__selectionType == "CASCADINGELITE"):
            selectionType = SelectionType.CASCADINGELITE
        elif (self.__selectionType == "ROULETTEELITE"):
            selectionType = SelectionType.ROULETTEELITE
        elif (self.__selectionType == "ROULETTETOTAL"):
            selectionType = SelectionType.ROULETTETOTAL
        else:
            self.__buildNotificationDialog("Selection type " + self.__selectionType + " is not properly defined.")
            return
        terminationType = -1
        if (self.__terminationType == "CONVERGENCE"):
            terminationType = TerminationType.CONVERGENCE
        elif (self.__terminationType == "OPTIMUMCONVERGENCE"):
            terminationType = TerminationType.OPTIMUMCONVERGENCE
        elif (self.__terminationType == "NATURAL"):
            terminationType = TerminationType.NATURAL
        else:
            self.__buildNotificationDialog("Termination type " + self.__terminationType + " is not properly defined.")
            return
        self.__algorithm = Algorithm(self.__numGenes, self.__popSize, self.__maxGen, crossoverType, self.__crossoverChance, self.__mutationSize,
                                     self.__mutationChance, self.__eliteSize, selectionType, terminationType, self.__optimumX,
                                     self.__optimumY, self.__convergenceConfidence, self.__evaluationDictionary)
        startTime = time.clock() * 1000
        bfRes = self.__algorithm.bruteForce()
        print("Milliseconds Elapsed:",(time.clock() * 1000) - startTime)
        print("Best Score:", bfRes[0])
        print("Locations:", bfRes[1])

    def startGUI(self):
        self.__root.mainloop()

    def destroyGUI(self):
        self.__root.destroy()


