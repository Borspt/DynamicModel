import json
from ElementsTechnological import *
from AdditionalFunctions import *
from matplotlib import pyplot as plt
import time
import os
import keyboard

JSON = "schemes/schemetest (3) (2pump).json"

with open(JSON, "r", encoding="utf-8") as load_file:
    test_json = json.load(load_file)

elements_dict = {}

PROCESSTIME = 1200
TAU = 1
SOUNDSPEED = 1000
DENSITY = 850
INITIAL_VELOCITY = 0.01
INITIAL_PRESSURE = 30 * g * 10000
SHOW_GRAPH = True
SAVE_LOGS = True
showSpeed = 100
USER_CHOICE = [123, 129]
graph_ylim = 3
startStep = 0
miniStepRatio = 100

all_objects, pipeObjects, boundaryObjects, branchIdList = init_objects(test_json=test_json)

steps = math.ceil(PROCESSTIME / TAU)

boundaryConditions = calcInitValues(boundaryObjects=boundaryObjects, initialVelocity=INITIAL_VELOCITY,
                                    initialPressure=INITIAL_PRESSURE, density=DENSITY)

for pipe in pipeObjects.values():
    pipe.initMesh(soundSpeed=SOUNDSPEED, steps=steps, tau=1, miniStepRatio=miniStepRatio)

    pipe.velocityMesh[0].fill(0.01)
    for index in range(len(pipe.pressureMesh[0])):
        pipe.pressureMesh[0][index] = INITIAL_PRESSURE - pipe.heights[index] * g * DENSITY
    if pipe.isTechnological:
        pipe.velocityMiniMesh[0].fill(0.01)
        pipe.pressureMiniMesh[0].fill(INITIAL_PRESSURE - pipe.heights[0] * g * DENSITY)

    boundaryLeftVelocity, boundaryLeftPressure, boundaryRightVelocity, boundaryRightPressure = \
        getBoundaryConditions(boundaryDict=boundaryConditions, boundaryObjects=boundaryObjects,
                              branchIdList=branchIdList, leftNeighbor=pipe.neighbors[0],
                              rightNeighbor=pipe.neighbors[1], pipeId=pipe.id)

    pipe.pressureMesh[0][0], pipe.pressureMesh[0][-1] = boundaryLeftPressure, boundaryRightPressure
    pipe.velocityMesh[0][0], pipe.velocityMesh[0][-1] = boundaryLeftVelocity, boundaryRightVelocity

    '''
    Заполнение
    граничных
    условий
    '''
calcTimeStart = time.time()
for stepNumber in range(1, steps + 1):
    # os.system('cls')
    # print(f'Шаг расчета {stepNumber} из {steps + 1}')
    # print('\n')
    # time.sleep(0.0001)

    for boundaryElement in boundaryObjects.values():
        if isinstance(boundaryElement, (FlowPressureSetter, Tank)):
            if boundaryElement.start:
                neighborId = boundaryElement.neighbors[1]
                neighborPoint = 1
            else:
                neighborId = boundaryElement.neighbors[0]
                neighborPoint = -2
            pipe = pipeObjects[neighborId]
            if pipe.isTechnological:
                pipeResistance = 0
            else:
                pipeResistance = pipe.calcResistance(velocity=(pipe.velocityMesh[stepNumber - 1][neighborPoint]))

            boundaryLeftPressure, boundaryRightPressure, \
                boundaryLeftVelocity, boundaryRightVelocity = \
                boundaryElement.calcBoundaries(tau=TAU, soundSpeed=SOUNDSPEED, pipeHeight=pipe.heights[neighborPoint],
                                               pipeVelocity=pipe.velocityMesh[stepNumber - 1][neighborPoint],
                                               pipePressure=pipe.pressureMesh[stepNumber - 1][neighborPoint],
                                               pipeResistance=pipeResistance,
                                               pipeKCorrection=pipe.kCorrection)
            conditions = {
                'pressure': [boundaryLeftPressure, boundaryRightPressure],
                'velocity': [boundaryLeftVelocity, boundaryRightVelocity]
            }
        elif isinstance(boundaryElement, Branch):
            pipeNeighbors = []
            pipeVelocityList = []
            pipePressureList = []
            pipeDiametersList = []
            pipeResistanceList = []
            pipeKCorrectionList = []
            pipeGravityFactorList = []

            for pipeId in boundaryElement.neighbors:
                pipe = pipeObjects[pipeId]
                pipeNeighbors.append(pipe)
                if pipeId in boundaryElement.neighborsOut:
                    pipeVelocity = pipe.velocityMesh[stepNumber - 1][1]
                    pipeVelocityList.append(pipeVelocity)
                    pipePressureList.append(pipe.pressureMesh[stepNumber - 1][1])
                    try:
                        pipeGravityFactorList.append((boundaryElement.height - pipe.start_height))
                    except Exception as e:
                        print(e)
                        pipeGravityFactorList.append(0)


                else:
                    pipeVelocity = pipe.velocityMesh[stepNumber - 1][-2]
                    pipeVelocityList.append(pipeVelocity)
                    pipePressureList.append(pipe.pressureMesh[stepNumber - 1][-2])
                    try:
                        pipeGravityFactorList.append((pipe.end_height - boundaryElement.height))
                    except Exception as e:
                        print(e)
                        pipeGravityFactorList.append(0)
                try:
                    if pipe.isTechnological:
                        pipeResistance = 0
                    else:
                        pipe.calcResistance(velocity=pipe.calcResistance(velocity=pipeVelocity))
                    pipeResistanceList.append(pipeResistance)
                except Exception as e:
                    print(e)
                pipeDiametersList.append(pipe.innerDiameter)
                pipeKCorrectionList.append(pipe.kCorrection)
            densityList = [DENSITY, DENSITY, DENSITY]
            soundSpeedList = [SOUNDSPEED, SOUNDSPEED, SOUNDSPEED]
            boundaryPressure, boundaryVelocity = \
                boundaryElement.calcBoundaries(tau=TAU, soundSpeed=soundSpeedList, density=densityList,
                                               gravityFactor=pipeGravityFactorList, velocity=pipeVelocityList,
                                               pressure=pipePressureList, diameter=pipeDiametersList,
                                               resistance=pipeResistanceList, kCorrection=pipeKCorrectionList,
                                               pipeObjects=pipeNeighbors)
            conditions = {
                'pressure': boundaryPressure,
                'velocity': boundaryVelocity
            }


        else:
            leftNeighborId = boundaryElement.neighbors[0]
            rightNeighborId = boundaryElement.neighbors[1]
            leftPipe = pipeObjects[leftNeighborId]
            rightPipe = pipeObjects[rightNeighborId]
            forwDiameter = rightPipe.innerDiameter
            prevDiameter = leftPipe.innerDiameter
            # if rightPipe.isTechnological:
            #     miniPrevVelocity = leftPipe.velocityMesh[stepNumber - 1][-2]
            #     miniPrevPressure = leftPipe.pressureMesh[stepNumber - 1][-2]
            #     # currConditions = boundaryConditions[boundaryElement.id]
            #     # boundaryLeftVelocity = currConditions['velocity'][1]
            #     # boundaryRightVelocity = currConditions['pressure'][1]
            #     # print(currConditions)
            #     for miniStep in range(1, miniStepRatio):
            #         miniForwVelocity = rightPipe.velocityMiniMesh[miniStep][1]
            #         miniForwPressure = rightPipe.pressureMiniMesh[miniStep][1]
            #         boundaryMiniLeftPressure, boundaryMiniRightPressure, boundaryMiniLeftVelocity, \
            #             boundaryMiniRightVelocity = \
            #             boundaryElement.calcBoundaries(tau=TAU, soundSpeed=SOUNDSPEED,
            #                                            density=DENSITY,
            #                                            forwG=0,
            #                                            prevG=0,
            #                                            forwVelocity=miniForwVelocity,
            #                                            prevVelocity=miniPrevVelocity,
            #                                            forwPressure=miniForwPressure,
            #                                            prevPressure=miniPrevPressure,
            #                                            forwPipeResistance=0,
            #                                            prevPipeResistance=leftPipe.calcResistance(
            #                                                velocity=miniPrevVelocity),
            #                                            forwKCorrection=0,
            #                                            prevKCorrection=leftPipe.kCorrection/miniStepRatio,
            #                                            prevDiameter=prevDiameter,
            #                                            forwDiameter=forwDiameter)
            #         pipe.velocityMiniMesh[miniStep]

            forwVelocity = rightPipe.velocityMesh[stepNumber - 1][1]
            prevVelocity = leftPipe.velocityMesh[stepNumber - 1][-2]
            forwPressure = rightPipe.pressureMesh[stepNumber - 1][1]
            prevPressure = leftPipe.pressureMesh[stepNumber - 1][-2]

            try:
                boundaryHeight = boundaryElement.height
                rightHeight = rightPipe.start_height
                leftHeight = leftPipe.end_height
                forwG = boundaryHeight - rightHeight
                prevG = leftHeight - boundaryHeight
            except Exception as e:
                print(e)
                forwG = 0
                prevG = 0

            boundaryLeftPressure, boundaryRightPressure, boundaryLeftVelocity, boundaryRightVelocity = \
                boundaryElement.calcBoundaries(tau=TAU, soundSpeed=SOUNDSPEED, density=DENSITY,
                                               forwG=forwG,
                                               prevG=prevG,
                                               forwVelocity=forwVelocity, prevVelocity=prevVelocity,
                                               forwPressure=forwPressure, prevPressure=prevPressure,
                                               forwPipeResistance=rightPipe.calcResistance(velocity=forwVelocity),
                                               prevPipeResistance=leftPipe.calcResistance(velocity=prevVelocity),
                                               forwKCorrection=rightPipe.kCorrection,
                                               prevKCorrection=leftPipe.kCorrection, prevDiameter=prevDiameter,
                                               forwDiameter=forwDiameter)

            conditions = {
                'pressure': [boundaryLeftPressure, boundaryRightPressure],
                'velocity': [boundaryLeftVelocity, boundaryRightVelocity]
            }

        boundaryConditions[boundaryElement.id] = conditions

    for pipe in pipeObjects.values():
        # print(pipe)
        leftNeighbor = pipe.neighbors[0]
        rightNeighbor = pipe.neighbors[1]
        boundaryLeftVelocity, boundaryLeftPressure, boundaryRightVelocity, boundaryRightPressure = \
            getBoundaryConditions(boundaryDict=boundaryConditions, boundaryObjects=boundaryObjects,
                                  branchIdList=branchIdList, leftNeighbor=leftNeighbor, rightNeighbor=rightNeighbor,
                                  pipeId=pipe.id)
        pipe.calcStep(soundSpeed=SOUNDSPEED, density=DENSITY, boundaryLeftVelocity=boundaryLeftVelocity,
                      boundaryRightVelocity=boundaryRightVelocity, boundaryLeftPressure=boundaryLeftPressure,
                      boundaryRightPressure=boundaryRightPressure, stepNumber=stepNumber, tau=TAU)
calcTimeEnd = time.time()
print(f'Расчет {steps + 1} шагов произведен за {calcTimeEnd - calcTimeStart} c')

showPauseTrigger = False


def revertPause():
    global showPauseTrigger
    if showPauseTrigger is False:
        showPauseTrigger = True
        print('Pause')
        time.sleep(1)
    else:
        showPauseTrigger = False
        print('Unpause')
        time.sleep(1)


def increaseStep():
    global step
    step += 1


def decreaseStep():
    global step
    step -= 1


if SAVE_LOGS:
    for elementId in USER_CHOICE:
        pipe = pipeObjects[elementId]
        velocity = pipe.velocityMesh
        np.savetxt(f"Logs/Pipe{elementId}.csv", velocity, delimiter=",")
        print(f'Логи по трубе {elementId} сохранены')

if SHOW_GRAPH:
    plt.ion()
    fig = plt.figure(1)
    figv = plt.figure(2)
    timeModel = 0
    x = []
    y = []
    x_end = []
    y_end = []
    for elementId in USER_CHOICE:
        pipe = pipeObjects[elementId]
        if pipe.isTechnological:
            pass
        else:
            for i in range(len(pipe.profile)):  ## Почему тут берутся данные по профилю, а не сетке?
                # - строится красный профиль высот и отмечаются концы труб

                el = pipe.profile[i]
                if i == 0 or (i == len(pipe.profile) - 1):
                    x_end.append(el['distance'])
                    y_end.append(el['height'])
                x.append(el['distance'])
                y.append(el['height'])
        # print(f'Y = {y}')
    step = startStep - 1
    keyboard.add_hotkey('space', revertPause)
    keyboard.add_hotkey('right', increaseStep)
    keyboard.add_hotkey('left', decreaseStep)
    while step < stepNumber:
        # event = keyboard.read_event()
        # keyboard.add_hotkey('space', lambda: print('space was pressed'))
        if showPauseTrigger == False:
            step += 1
        else:
            pass

        timeModel = TAU * step
        dotX = np.array([])
        ph = np.array([])
        velocity = np.array([])
        heights = np.array([])
        for elementId in USER_CHOICE:
            pipe = pipeObjects[elementId]
            additionalX = np.array(pipe.meshX) / 1000 + pipe.profile[0]['distance']
            dotX = np.append(dotX, additionalX)
            heights = np.append(heights, pipe.heights)

            for i in range(len(pipe.pressureMesh[0])):
                # ph.append((self.pressureMesh[step-1][i] + density*g*heights[i])/density/g)
                additionalPh = (pipe.pressureMesh[step - 1][i] + DENSITY * g * pipe.heights[i]) / DENSITY / g
                ph = np.append(ph, additionalPh)
            for i in range(len(pipe.velocityMesh[0])):
                # ph.append((self.pressureMesh[step-1][i] + density*g*heights[i])/density/g)
                velocity = np.append(velocity, pipe.velocityMesh[step - 1][i])
        # print(ph)
        time_start = time.time()
        plt.figure(1)
        plt.clf()
        plt.plot(dotX, ph, 'blue', marker='o', label='Гидроуклон')
        plt.plot(x_end, y_end, 'black', marker='D', label='Концы')
        plt.plot(x, y, color='red', ls='--', label='Высота пролегания')
        plt.scatter(dotX, heights, marker='o')
        step_text = f'Шаг: {step} из {stepNumber + 1}'
        plt.text(0.95, 0.95, step_text, transform=plt.gca().transAxes, fontsize=10,
                 verticalalignment='top', horizontalalignment='right')
        plt.title(f'Обновляемый график давления {timeModel:.1f} с')
        plt.xlabel('Дистанция, км', fontsize=30)
        plt.ylabel('Напор, м', fontsize=30)
        plt.yticks(fontsize=20)
        plt.xticks(fontsize=20)

        plt.figure(2)
        plt.clf()
        plt.plot(dotX, velocity, 'blue', marker='o', label='Скорость')

        # plt.scatter(dotX, self.heights, marker='o')
        step_text = f'Шаг: {step} из {stepNumber + 1}'
        plt.text(0.95, 0.95, step_text, transform=plt.gca().transAxes, fontsize=10,
                 verticalalignment='top', horizontalalignment='right')
        plt.title(f'Обновляемый график скорости {timeModel:.1f} с')
        plt.xlabel('Дистанция, км', fontsize=30)
        plt.ylabel('Скорость , м/с', fontsize=30)
        plt.ylim(top=graph_ylim, bottom=-graph_ylim)
        plt.yticks(fontsize=20)
        plt.xticks(fontsize=20)
        end_time = time.time()
        graph_time = end_time - time_start
        # print(f'Время на отображение {graph_time} c')
        if (TAU / showSpeed - graph_time) > 0:
            plt.pause(TAU / showSpeed - graph_time)
        else:
            plt.pause(0.00001)
        # plt.pause(1)

    # plt.savefig('graph (3 tank + 2 fps) (20Rotor).png')

#
# print(pipe.__dict__)
# print(pipe.profile)
# print(pipe.heights)
# print(pipe.meshX)
# print(len(pipe.heights))
# print(len(pipe.meshX))
#
# print('_________________________________________________________________')
# print('Velocity Mesh:')
# print(pipe.velocityMesh)
# print('_________________________________________________________________')
# print(pipe.pressureMesh)
