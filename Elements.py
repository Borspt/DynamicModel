import csv
import math
import numpy as np
import warnings
import os

g = 9.81


def invariantMinus(forwPressure, density, soundSpeed, forwVelocity, kCorrection, tau, forwResistance, gravityFactor):
    invariant = forwPressure - density * soundSpeed * forwVelocity + density * g * gravityFactor \
                + soundSpeed * tau * kCorrection * (density * forwResistance)
    return invariant


def invariantPlus(prevPressure, density, soundSpeed, prevVelocity, kCorrection, tau, prevResistance, gravityFactor):
    invariant = prevPressure + density * soundSpeed * prevVelocity - density * g * gravityFactor \
                - soundSpeed * tau * kCorrection * (density * prevResistance)
    return invariant


def getBoundaryConditions(boundaryDict, boundaryObjects, leftNeighbor, rightNeighbor, pipeId):
    boundaryLeft = boundaryDict.get(leftNeighbor, None)
    boundaryRight = boundaryDict.get(rightNeighbor, None)

    assert boundaryLeft is not None, f'Pipe id{pipeId} left neighbor not in boundaryConditions'
    assert boundaryRight is not None, f'Pipe id{pipeId} right neighbor not in boundaryConditions'

    if isinstance(boundaryObjects[leftNeighbor], Branch):
        boundaryLeftVelocity = None
        boundaryLeftPressure = None
        branch = boundaryObjects[leftNeighbor]
        for index in range(len(branch.neighbors)):
            if pipeId == branch.neighbors[index]:
                boundaryLeftVelocity = boundaryLeft.get('velocity', None)[index]
                boundaryLeftPressure = boundaryLeft.get('pressure', None)[index]
            else:
                pass

    else:
        boundaryLeftVelocity = boundaryLeft.get('velocity', None)[1]
        boundaryLeftPressure = boundaryLeft.get('pressure', None)[1]

    assert boundaryLeftVelocity is not None, f'Pipe id{pipeId} left velocity is None, leftNeighbor id{leftNeighbor}'
    assert boundaryLeftPressure is not None, f'Pipe id{pipeId} left pressure is None, leftNeighbor id{leftNeighbor}'

    if isinstance(boundaryObjects[rightNeighbor], Branch):
        boundaryRightVelocity = None
        boundaryRightPressure = None
        branch = boundaryObjects[rightNeighbor]
        for index in range(len(branch.neighbors)):
            if pipeId == branch.neighbors[index]:
                boundaryRightVelocity = boundaryRight.get('velocity', None)[index]
                boundaryRightPressure = boundaryRight.get('pressure', None)[index]
            else:
                pass
    else:
        boundaryRightVelocity = boundaryRight.get('velocity', None)[0]
        boundaryRightPressure = boundaryRight.get('pressure', None)[0]

    assert boundaryRightVelocity is not None, f'Pipe id{pipeId} right velocity is None, rightNeighbor id{rightNeighbor}'
    assert boundaryRightPressure is not None, f'Pipe id{pipeId} right pressure is None, rightNeighbor id{rightNeighbor}'

    return boundaryLeftVelocity, boundaryLeftPressure, boundaryRightVelocity, boundaryRightPressure


# J_minus = forwPressure - density * soundSpeed * forwVelocity +\
#                        density * g * forwG +\
#                       soundSpeed * tau * self.kCorrection * (density *forwB)
# J_plus = prevPressure + density * soundSpeed * prevVelocity\
#                      - density * g * prevG -\
#                      soundSpeed * tau * self.kCorrection * (density *prevB)


class BlockValve:
    def __init__(self, elementId, neighbors, height, resistanceCoefficient=0, status=1, fullFlow=1, kvFile=None):
        self.id = elementId
        self.fullFlow = fullFlow
        self.status = status
        if kvFile is not None:
            self.UploadKvFile(kvFile)
        self.kvFile = kvFile
        self.neighbors = neighbors
        self.height = height
        self.resistanceCoefficient = resistanceCoefficient
    def deltaP(self, density, flow):
        self.deltaPressure = 0
        if self.status != 'open':
            raise ValueError('Incorrect Schema. DeltaP at closed blockValve')
        else:
            deltaP = 0
            return deltaP

    def UploadKvFile(self, kvFile=None):
        if kvFile is not None:
            self.kvFile = kvFile
        if self.kvFile is not None:
            with open(self.kvFile, 'r', newline='') as file:
                kvCharacteristic = csv.reader(file, delimiter=';')
                print(kvCharacteristic)
                kvChars = {}
                for row in kvCharacteristic:
                    try:
                        # print(type(row[0]))
                        kvChars[float(row[0])] = float(row[1])
                    except Exception as e:
                        print(e)
                    if isinstance(row[0], float):
                        print(row[0])
                print(kvChars)
                self.Chars = kvChars

    def calcFlow(self, openPercentage):
        for index in range(len(self.Chars.keys()) - 1):
            if (openPercentage >= list(self.Chars.keys())[index]) and (
                    openPercentage <= list(self.Chars.keys())[index + 1]):
                key1 = list(self.Chars.keys())[index]
                key2 = list(self.Chars.keys())[index + 1]
                flowPercentage = self.Chars[key1] * (key2 - openPercentage) / (key2 - key1) \
                                 + self.Chars[key2] * (openPercentage - key1) / (key2 - key1)
                curFlow = self.fullFlow * flowPercentage
                print(curFlow)
                return curFlow

    def calcBoundaries(self, tau, soundSpeed, density, forwG, prevG, forwVelocity, prevVelocity, forwPressure,
                       prevPressure, prevPipeResistance, forwPipeResistance, forwKCorrection, prevKCorrection,
                       prevDiameter, forwDiameter):

        Jminus = invariantMinus(forwPressure=forwPressure, density=density, soundSpeed=soundSpeed,
                                forwVelocity=forwVelocity, kCorrection=forwKCorrection, tau=tau,
                                forwResistance=forwPipeResistance, gravityFactor=forwG)
        Jplus = invariantPlus(prevPressure=prevPressure, density=density, soundSpeed=soundSpeed,
                              prevVelocity=prevVelocity, kCorrection=prevKCorrection, tau=tau,
                              prevResistance=prevPipeResistance, gravityFactor=prevG)

        velocity = (Jplus - Jminus) / (density * soundSpeed) / (
                1 + math.sqrt(1 + self.resistanceCoefficient / 2 * density * soundSpeed * (Jplus - Jminus)))
        boundaryLeftVelocity = velocity
        boundaryRightVelocity = velocity

        boundaryLeftPressure = boundaryRightPressure = (Jplus + Jminus) / 2

        return boundaryLeftPressure, boundaryRightPressure, boundaryLeftVelocity, boundaryRightVelocity


class FilterStrainer:
    def __init__(self, elementiId, neighbors, flow=0):
        self.id = elementiId
        self.flow = flow
        self.userChoice = 'noLosses'
        self.deltaPressure = 0
        self.neighbors = neighbors

    def deltaP(self, density, flow):
        if self.userChoice == 'noLosses':
            deltaP = 0
            return deltaP
        else:
            print('No method for filterStrainer // TBD')
            deltaP = 0
            return deltaP


class Pipe:
    def __init__(self, elementId, innerDiameter, length, profile, start_height, end_height, neighbors, technoNeighbors,
                 isTechnological, technologicalResistance=0, roughness=0.0002, losses='auto', viscosity=1e-3):

        self.velocityMesh = None
        self.id = elementId
        self.innerDiameter = innerDiameter
        self.length = length * 1000
        self.start_height = start_height
        self.end_height = end_height
        self.viscosity = viscosity
        self.rheinolds = None
        self.losses = losses
        self.roughness = roughness
        self.profile = profile
        self.neighbors = neighbors
        self.isTechnological = isTechnological
        self.technologicalResistance = technologicalResistance
        self.technoNeighbors = technoNeighbors

    def calcLambda(self, flow, useFlow=True):
        if useFlow:
            self.velocity = 4 * abs(flow) / (math.pi * pow(self.innerDiameter, 2))
        else:
            self.velocity = flow
        try:
            self.rheinolds = abs(self.velocity) * self.innerDiameter / self.viscosity
        except:
            print('Id', self.id)
        # print(self.rheinolds)
        if self.rheinolds == 0:
            print('Rheinolds = 0 Id', self.id)
            self.rheinolds = 0.001

        if (self.losses == 'Stokes') or ((self.losses == 'auto') and (self.rheinolds <= 2000)):
            try:
                _lambda = 64 / self.rheinolds
            except ZeroDivisionError:
                _lambda = 64 / (self.rheinolds + 0.00000000000000001)
        if (self.losses == 'Blasius') or (
                (self.losses == 'auto') and (self.rheinolds >= 2000) and (self.rheinolds <= 10 / self.roughness)):
            _lambda = 0.3164 / pow(self.rheinolds, 0.25)
        if (self.losses == 'Altshul') or ((self.losses == 'auto') and (self.rheinolds >= 10 / self.roughness) and (
                self.rheinolds <= 500 / self.roughness)):
            _lambda = 0.11 * pow(68 / self.rheinolds + self.roughness, 0.25)
        if (self.losses == 'Isaev') or ((self.losses == 'auto') and (self.rheinolds >= 10 / self.roughness) and (
                self.rheinolds <= 500 / self.roughness)):
            _lambda = 1 / (1.18 * math.log(6.8 / self.rheinolds + pow(self.roughness / 3.7, 1.11))) ** 2
        if (self.losses == 'Ginzburg') or (
                (self.losses == 'auto') and (self.rheinolds >= 2000) and (self.rheinolds <= 10000)):
            y = 1 - math.exp(-0.002 * (self.rheinolds - 2320))
            _lambda = 0.3164 / pow(self.rheinolds, 0.25) * y + 64 * (1 - y) / self.rheinolds
        if (self.losses == 'Nikuradze') or ((self.losses == 'auto') and (self.rheinolds >= 500 / self.roughness)):
            _lambda = 1 / (1.14 - 21 * g * self.roughness) ** 2
        if (self.losses == 'Shifrinson') or ((self.losses == 'auto') and (self.rheinolds >= 500 / self.roughness)):
            _lambda = 0.11 * pow(self.roughness, 0.25)
        try:
            return _lambda
        except:
            print(self.rheinolds)

    def calcResistance(self, velocity):
        if self.isTechnological:
            resistance = self.technologicalResistance
        else:
            _lambda = self.calcLambda(velocity, useFlow=False)
            try:
                resistance = _lambda * velocity * abs(velocity) / (2 * self.innerDiameter)
            except Exception as e:
                print(e)
        return resistance

    def deltaP(self, density, flow):
        global g
        _lambda = self.calcLambda(flow)
        height_pressure = g * (self.end_height - self.start_height) * density
        self.density = density
        self.height_pressure = height_pressure
        # self.losses = -_lambda * self.velocity * abs(self.velocity) * self.length * density / (2 * self.innerDiameter)
        deltaP = -_lambda * self.velocity * abs(self.velocity) * self.length * density / (2 * self.innerDiameter) \
                 - height_pressure
        self.deltaPressure = deltaP
        return deltaP

    def __calcMeshX__(self, soundSpeed, tau):
        if self.length > 0.5 * soundSpeed * tau:
            self.meshCountX = math.floor(self.length / soundSpeed / tau) + 1
        else:
            self.meshCountX = 2

        self.lengthCalc = (self.meshCountX - 1) * soundSpeed * tau
        self.kCorrection = self.length / self.lengthCalc
        self.meshDeltaX = soundSpeed * tau * self.kCorrection
        self.meshX = [x * self.meshDeltaX for x in range(self.meshCountX - 1)]
        self.meshX.append(self.length)

    def __interpolation__(self, xPoint):
        profile = self.profile
        i = 0
        xPoint = xPoint / 1000
        while i != len(profile):
            x = profile[i]['distance']
            x0 = profile[0]['distance']

            if x > xPoint + x0:
                a = (profile[i]['height'] - profile[i - 1]['height']) / (
                        profile[i]['distance'] - profile[i - 1]['distance'])
                b = profile[i]['height'] - a * profile[i]['distance']
                height = a * (xPoint + x0) + b
                return height
            i += 1
        return None

    def initMesh(self, soundSpeed, steps, tau):
        '''
        Определение координат X для расчета
        '''
        self.__calcMeshX__(soundSpeed, tau=tau)
        '''
        Расчет высот соответствующих координате X
        '''
        heights = []
        # heights = [321.95] * len(self.meshX)
        if self.isTechnological:
            heights.append(self.start_height)
            heights.append(self.end_height)
            self.heights = heights
        else:
            heights.append(self.profile[0]['height'])
            for xPoint in self.meshX[1:-1]:
                height = self.__interpolation__(xPoint)
                if height is None:
                    print(f'None height pipeId{self.id}')
                heights.append(height)
            heights.append(self.profile[-1]['height'])
            self.heights = heights

        '''
        Инициализация сетки
        '''
        self.velocityMesh = np.zeros((steps + 1, len(self.heights)))
        self.pressureMesh = np.zeros((steps + 1, len(self.heights)))

        return

    def calcStep(self, soundSpeed, density, boundaryLeftPressure, boundaryRightPressure,
                 boundaryLeftVelocity, boundaryRightVelocity, stepNumber, tau):

        '''
        Расчет для случая с резервуаром, устарел
        '''

        # self.velocityMesh[stepNumber][0] = 1 / (density * soundSpeed) * (
        #         boundaryLeft + density * soundSpeed * self.velocityMesh[stepNumber - 1][1]
        #         - self.pressureMesh[stepNumber - 1][1] - tau * soundSpeed * self.kCorrection * (
        #                     density * g * (self.heights[0] - self.heights[1]) / soundSpeed * tau
        #                     + density * self.calcResistance(velocity=(self.velocityMesh[stepNumber - 1][1]))))
        #
        # self.velocityMesh[stepNumber][-1] = 1 / (density * soundSpeed) * (self.pressureMesh[stepNumber - 1][-2] +
        #
        #                                                             density * soundSpeed *
        #                                                             self.velocityMesh[stepNumber - 1][
        #                                                                 -2] - boundaryRight -
        #                                                             soundSpeed * tau * self.kCorrection * (
        #                                                                         density * g * (
        #                                                                             self.heights[-1] - self.heights[
        #                                                                         -2]) / soundSpeed * tau
        #                                                                         + density * self.calcResistance(
        #                                                                     velocity=(
        #                                                                     self.velocityMesh[stepNumber - 1][-2]))))
        '''
                Заполнение краевых условий
                '''
        self.velocityMesh[stepNumber][0] = boundaryLeftVelocity
        self.velocityMesh[stepNumber][-1] = boundaryRightVelocity
        self.pressureMesh[stepNumber][0] = boundaryLeftPressure
        self.pressureMesh[stepNumber][-1] = boundaryRightPressure

        for i in range(1, self.pressureMesh.shape[1] - 1):
            prevPressure = self.pressureMesh[stepNumber - 1][i - 1]
            forwPressure = self.pressureMesh[stepNumber - 1][i + 1]
            prevVelocity = self.velocityMesh[stepNumber - 1][i - 1]
            curVelocity = self.velocityMesh[stepNumber - 1][i]
            forwVelocity = self.velocityMesh[stepNumber - 1][i + 1]

            forwG = self.heights[i + 1] - self.heights[i]
            prevG = self.heights[i] - self.heights[i - 1]
            if self.isTechnological:
                forwB = 0
                prevB = 0
            else:
                forwB = self.calcResistance(velocity=(forwVelocity))
                prevB = self.calcResistance(velocity=(prevVelocity))

            if self.isTechnological:

                Jminus = forwPressure - density * soundSpeed * forwVelocity + density * g * forwG

                Jplus = prevPressure + density * soundSpeed * prevVelocity - density * g * prevG


            else:

                Jminus = invariantMinus(forwPressure=forwPressure, density=density, soundSpeed=soundSpeed,
                                        forwVelocity=forwVelocity, kCorrection=self.kCorrection, tau=tau,
                                        forwResistance=forwB, gravityFactor=forwG)
                Jplus = invariantPlus(prevPressure=prevPressure, density=density, soundSpeed=soundSpeed,
                                      prevVelocity=prevVelocity, kCorrection=self.kCorrection, tau=tau,
                                      prevResistance=prevB, gravityFactor=prevG)

            self.pressureMesh[stepNumber][i] = (Jplus + Jminus) / 2
            self.velocityMesh[stepNumber][i] = (Jplus - Jminus) / (2 * density * soundSpeed)

            # prevB = self.calcResistance(velocity=prevVelocity)
            # forwB = self.calcResistance(velocity=(forwVelocity))
            # prevG = g * (heights[i] - heights[i - 1])
            # forwG = g * (heights[i+1] - heights[i])
            # self.velocityMesh[step][i] = 1 / (soundSpeed * 2 * density) * (prevPressure - forwPressure +
            #                                                         density * soundSpeed * prevVelocity + density * soundSpeed * forwVelocity
            #                                                         - soundSpeed * tau * self.kCorrection * (density * forwB + density * prevB) -
            #                                                         (density * prevG + density * forwG))
            #
            # self.pressureMesh[step][i] = density * soundSpeed * self.velocityMesh[step][i] + forwPressure - \
            #                           density * soundSpeed * forwVelocity \
            #                           + density * soundSpeed * tau*self.kCorrection*prevB + density * forwG


class CheckValve:
    def __init__(self, elementId, neighbors, height, resistanceCoefficient=0, flow=0):
        self.id = elementId
        self.flow = flow
        self.deltaPressure = 0
        self.height = height
        self.neighbors = neighbors
        self.resistanceCoefficient = resistanceCoefficient

    def deltaP(self, density, flow):

        self.flow = flow

        if flow >= -0:
            deltaP = 0
        else:
            deltaP = 10000000 * (math.exp(-self.flow / 1.39) - 1)
        self.deltaPressure = deltaP
        return deltaP

    def calcBoundaries(self, tau, soundSpeed, density, forwG, prevG, forwVelocity, prevVelocity, forwPressure,
                       prevPressure, prevPipeResistance, forwPipeResistance, forwKCorrection, prevKCorrection,
                       prevDiameter, forwDiameter):

        Jminus = invariantMinus(forwPressure=forwPressure, density=density, soundSpeed=soundSpeed,
                                forwVelocity=forwVelocity, kCorrection=forwKCorrection, tau=tau,
                                forwResistance=forwPipeResistance, gravityFactor=forwG)
        Jplus = invariantPlus(prevPressure=prevPressure, density=density, soundSpeed=soundSpeed,
                              prevVelocity=prevVelocity, kCorrection=prevKCorrection, tau=tau,
                              prevResistance=prevPipeResistance, gravityFactor=prevG)

        velocity = (Jplus - Jminus) / (density * soundSpeed) / (
                1 + math.sqrt(1 + self.resistanceCoefficient / 2 * density * soundSpeed * (Jplus - Jminus)))
        if velocity >= 0:
            boundaryLeftVelocity = velocity
            boundaryRightVelocity = velocity
            boundaryLeftPressure = boundaryRightPressure = (Jplus + Jminus) / 2
        else:
            boundaryLeftVelocity = 0
            boundaryRightVelocity = 0
            boundaryLeftPressure = Jplus
            boundaryRightPressure = Jminus



        return boundaryLeftPressure, boundaryRightPressure, boundaryLeftVelocity, boundaryRightVelocity


class Pump:
    rotor_objects = None

    def __init__(self, elementId, neighbors, useConstHead, deltaHead, rotorId, height=None):
        self.id = elementId
        self.rotorId = rotorId
        self.useConstHead = useConstHead
        self.deltaHead = deltaHead
        self.neighbors = neighbors
        self.height = height

    def find_rotor(self):

        for rotor in self.rotor_objects:
            if rotor.id == self.rotorId:
                return rotor
        return None

    def deltaP(self, density, flow):
        rotor = self.find_rotor()
        deltaP = rotor.deltaP(density, flow)
        self.flow = flow
        self.deltaPressure = deltaP
        return deltaP

    def calcBoundaries(self, tau, soundSpeed, density, forwG, prevG, forwVelocity, prevVelocity, forwPressure,
                       prevPressure,
                       forwPipeResistance, prevPipeResistance, forwKCorrection, prevKCorrection, forwDiameter,
                       prevDiameter):
        prevArea = math.pi * prevDiameter ** 2 / 4
        forwArea = math.pi * forwDiameter ** 2 / 4

        Jminus = invariantMinus(forwPressure=forwPressure, density=density, soundSpeed=soundSpeed,
                                forwVelocity=forwVelocity, kCorrection=forwKCorrection, tau=tau,
                                forwResistance=forwPipeResistance, gravityFactor=forwG)
        Jplus = invariantPlus(prevPressure=prevPressure, density=density, soundSpeed=soundSpeed,
                              prevVelocity=prevVelocity, kCorrection=prevKCorrection, tau=tau,
                              prevResistance=prevPipeResistance, gravityFactor=prevG)

        if self.useConstHead == True:

            boundaryLeftVelocity = (self.deltaHead * density * g - (Jminus - Jplus)) / 2 / density / soundSpeed
            boundaryRightVelocity = boundaryLeftVelocity * prevArea / forwArea

            boundaryLeftPressure = Jplus - density * soundSpeed * boundaryLeftVelocity
            boundaryRightPressure = Jminus + density * soundSpeed * boundaryRightVelocity

            return boundaryLeftPressure, boundaryRightPressure, boundaryLeftVelocity, boundaryRightVelocity

        else:
            warnings.warn(message='Not const head // TBD')


class Rotor:
    def __init__(self, elementId, nominal_frequency, real_frequency, user_choice, density=None, flow=None):
        self.id = elementId
        self.user_choice = user_choice
        self.nominal_frequency = nominal_frequency
        self.real_frequency = real_frequency
        self.density = density
        self.flow = flow

    def deltaP(self, density, flow):
        global g
        deltaP = 0
        power = 0
        for coef in self.user_choice:
            deltaP += pow(self.real_frequency / self.nominal_frequency,
                          2) * g * density * coef * pow(flow * self.nominal_frequency / self.real_frequency, power)
            power += 1
        return deltaP


class Branch:
    def __init__(self, elementId, height, neighborsIn, neighborsOut):
        self.id = elementId
        self.deltaPressure = 0
        self.height = height
        self.neighborsIn = neighborsIn
        self.neighborsOut = neighborsOut
        self.neighbors = neighborsIn + neighborsOut

    def deltaP(self, density, flow):
        return 0

    def calcBoundaries(self, tau, soundSpeed, density, gravityFactor, velocity, pressure, diameter,
                       resistance, kCorrection, pipeObjects):
        pipeOutMask = [0, 0, 0]

        calcTwoPipeFlow = False
        pipeFlowIndex = []
        boundaryPressure = [None, None, None]
        boundaryVelocity = [None, None, None]
        directions = [None, None, None]
        # Шаг 3
        characteristic = [None, None, None]
        for i in range(len(self.neighbors)):
            pipe = pipeObjects[i]
            if pipe.id in self.neighborsOut:
                directions[i] = -1
                characteristic[i] = invariantMinus(forwPressure=pressure[i], density=density[i],
                                                   soundSpeed=soundSpeed[i],
                                                   forwVelocity=velocity[i], kCorrection=kCorrection[i], tau=tau,
                                                   forwResistance=resistance[i], gravityFactor=gravityFactor[i])
            else:
                directions[i] = 1
                characteristic[i] = invariantPlus(prevPressure=pressure[i], density=density[i],
                                                  soundSpeed=soundSpeed[i],
                                                  prevVelocity=velocity[i], kCorrection=kCorrection[i], tau=tau,
                                                  prevResistance=resistance[i], gravityFactor=gravityFactor[i])

        boundaryVelocity[2] = directions[2] * ((density[1] * soundSpeed[1]
                                                * math.pi * diameter[0] ** 2 / 4 * (
                                                        characteristic[2] - characteristic[0]) +
                                                density[0] * soundSpeed[0] * math.pi * diameter[1] ** 2 / 4 * (
                                                        characteristic[2] - characteristic[1])) /
                                               (density[0] * soundSpeed[0] * density[1] * soundSpeed[1] * math.pi *
                                                diameter[2] ** 2 / 4 +
                                                density[1] * soundSpeed[1] * density[2] * soundSpeed[2] * math.pi *
                                                diameter[0] ** 2 / 4 +
                                                density[0] * soundSpeed[0] * density[2] * soundSpeed[2] * math.pi *
                                                diameter[1] ** 2 / 4))

        boundaryPressure[0] = boundaryPressure[1] = boundaryPressure[2] = characteristic[2] - \
                                                                          directions[2] * density[2] * soundSpeed[2] * \
                                                                          boundaryVelocity[2]

        boundaryVelocity[0] = (characteristic[0] - boundaryPressure[0]) / (directions[0] * density[0] * soundSpeed[0])

        boundaryVelocity[1] = (characteristic[1] - boundaryPressure[1]) / (directions[1] * density[1] * soundSpeed[1])

        return boundaryPressure, boundaryVelocity


class Tank:
    def __init__(self, elementId, height, neighbors, inHeight=0, start=True, density=None, flow=0):
        self.id = elementId
        self.density = density
        self.flow = flow
        self.start = start
        self.inHeight = inHeight
        self.height = height
        self.neighbors = neighbors

    def deltaP(self, density=None):
        if density is None:
            density = self.density
        global g
        try:
            deltaP = (self.inHeight) * density * g
        except:
            print(self.id)

        # return deltaP
        if self.start:
            self.deltaPressure = deltaP
            return deltaP
        else:
            self.deltaPressure = -deltaP
            return -deltaP

    def calcBoundaries(self, tau, soundSpeed, pipeHeight, pipeVelocity, pipePressure, pipeResistance, pipeKCorrection):

        if self.start:
            boundaryLeftPressure = None
            boundaryRightPressure = self.deltaP()

            boundaryLeftVelocity = None
            Jminus = invariantMinus(forwPressure=pipePressure, density=self.density, soundSpeed=soundSpeed
                                    , forwVelocity=pipeVelocity, kCorrection=pipeKCorrection, tau=tau,
                                    forwResistance=pipeResistance, gravityFactor=(pipeHeight - self.height))

            boundaryRightVelocity = (boundaryRightPressure - Jminus) / (self.density * soundSpeed)


        else:

            boundaryLeftPressure = -self.deltaP()
            boundaryRightPressure = None

            Jplus = invariantPlus(prevPressure=pipePressure, density=self.density, soundSpeed=soundSpeed
                                  , prevVelocity=pipeVelocity, kCorrection=pipeKCorrection, tau=tau,
                                  prevResistance=pipeResistance, gravityFactor=(self.height - pipeHeight))

            boundaryLeftVelocity = (Jplus - boundaryLeftPressure) / (self.density * soundSpeed)

            boundaryRightVelocity = None

        return boundaryLeftPressure, boundaryRightPressure, boundaryLeftVelocity, boundaryRightVelocity


class FlowPressureSetter:
    def __init__(self, elementId, height, density, neighbors, start=True, dP=0, inHeight=0, flow=0):
        self.id = elementId
        self.density = density
        self.flow = flow
        self.inHeight = inHeight
        self.deltaPressure = dP
        self.height = height
        self.neighbors = neighbors
        self.start = start

    def deltaP(self, density=None):
        if density is None:
            density = self.density
        global g
        try:
            deltaP = (self.inHeight) * density * g
        except:
            print(self.id)

        # return deltaP
        if self.start:
            self.deltaPressure = deltaP
            return deltaP
        else:
            self.deltaPressure = -deltaP
            return -deltaP

    def calcBoundaries(self, tau, soundSpeed, pipeHeight, pipeVelocity, pipePressure, pipeResistance, pipeKCorrection):

        if self.start:
            boundaryLeftPressure = None
            boundaryRightPressure = self.deltaP()

            boundaryLeftVelocity = None
            boundaryRightVelocity = 1 / (self.density * soundSpeed) * (
                    boundaryRightPressure + self.density * soundSpeed * pipeVelocity
                    - pipePressure - tau * soundSpeed * pipeKCorrection * (
                            self.density * g * (self.height - pipeHeight) / soundSpeed / tau
                            + self.density * pipeResistance))
        else:

            boundaryLeftPressure = -self.deltaP()
            boundaryRightPressure = None

            boundaryLeftVelocity = 1 / (self.density * soundSpeed) * (pipePressure
                                                                      + self.density * soundSpeed * pipeVelocity - boundaryLeftPressure
                                                                      - soundSpeed * tau * pipeKCorrection * (
                                                                              self.density * g * (self.height
                                                                                                  - pipeHeight) / soundSpeed / tau + self.density * pipeResistance))
            boundaryRightVelocity = None
        return boundaryLeftPressure, boundaryRightPressure, boundaryLeftVelocity, boundaryRightVelocity


class Regulator:
    def __init__(self, elementId, neighbors, KP, KI, KD, target=0):
        self.id = elementId
        self.kp = KP
        self.ki = KI
        self.kd = KD
        self.sp = target
        self.error_last = 0
        self.integral_error = 0
        self.saturation_max = None
        self.saturation_min = None
        self.neighbors = neighbors

    def compute(self, pressure, dt):
        error = self.sp - pressure  # compute the error
        derivative_error = (
                                   error - self.error_last) / dt  # find the derivative of the error (how the error changes with time)
        self.integral_error += error * dt  # error build up over time
        output = self.kp * error + self.ki * self.integral_error + self.kd * derivative_error
        self.error_last = error
        if output > self.saturation_max and self.saturation_max is not None:
            output = self.saturation_max
        elif output < self.saturation_min and self.saturation_min is not None:
            output = self.saturation_min
        return output

    def setLims(self, min, max):
        self.saturation_max = max
        self.saturation_min = min


def interPolateVelocity(Vn, Vk, miniSteps, miniStep):
    return Vk + (Vn - Vk) * (miniSteps - miniStep) / miniSteps


def interPolatePressure(Pn, Pk, miniSteps, miniStep):
    return Pk + (Pn - Pk) * (miniSteps - miniStep) / miniSteps


class TechnologicalObject:
    def __init__(self, elementId, boundaryObjects, pipeObjects, miniSteps, tau, height, boundaryConditions={}):
        self.id = elementId
        self.boundaryObjects = boundaryObjects
        self.pipeObjects = pipeObjects
        self.miniSteps = miniSteps
        self.tau = tau
        self.boundaryConditions = boundaryConditions
        self.height = height

    def calcBoundaries(self, soundSpeed, density, stepNumber, globalPipeObjects):
        boundaryLeftP = None
        boundaryLeftV = None
        boundaryRightP = None
        boundaryRightV = None

        for miniStep in range(1, self.miniSteps + 1):
            for boundaryElement in self.boundaryObjects.values():
                if isinstance(boundaryElement, (FlowPressureSetter, Tank)):
                    if boundaryElement.start:
                        neighborId = boundaryElement.neighbors[1]
                        neighborPoint = 1
                    else:
                        neighborId = boundaryElement.neighbors[0]
                        neighborPoint = -2

                    pipe = self.pipeObjects[neighborId]

                    boundaryLeftPressure, boundaryRightPressure, \
                        boundaryLeftVelocity, boundaryRightVelocity = \
                        boundaryElement.calcBoundaries(tau=self.tau, soundSpeed=soundSpeed,
                                                       pipeHeight=pipe.heights[neighborPoint],
                                                       pipeVelocity=pipe.velocityMesh[stepNumber - 1][neighborPoint],
                                                       pipePressure=pipe.pressureMesh[stepNumber - 1][neighborPoint],
                                                       pipeResistance=0,
                                                       pipeKCorrection=0)
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
                        pipe = self.pipeObjects[pipeId]
                        pipeNeighbors.append(pipe)
                        if pipeId in boundaryElement.neighborsOut:
                            pipeVelocity = pipe.velocityMesh[stepNumber - 1][1]
                            pipeVelocityList.append(pipeVelocity)
                            pipePressureList.append(pipe.pressureMesh[stepNumber - 1][1])
                            try:
                                pipeGravityFactorList.append((pipe.heights[1] - boundaryElement.height))
                            except Exception as e:
                                print(e)
                                pipeGravityFactorList.append(0)


                        else:
                            pipeVelocity = pipe.velocityMesh[stepNumber - 1][-2]
                            pipeVelocityList.append(pipeVelocity)
                            pipePressureList.append(pipe.pressureMesh[stepNumber - 1][-2])
                            try:
                                pipeGravityFactorList.append((boundaryElement.height - pipe.heights[-2]))
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
                    densityList = [density, density, density]
                    soundSpeedList = [soundSpeed, soundSpeed, soundSpeed]
                    boundaryPressure, boundaryVelocity = \
                        boundaryElement.calcBoundaries(tau=self.tau, soundSpeed=soundSpeedList, density=densityList,
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
                    leftPipe = self.pipeObjects.get(leftNeighborId, None)
                    rightPipe = self.pipeObjects.get(rightNeighborId, None)
                    if leftNeighborId not in list(self.pipeObjects.keys()):
                        leftGlobalPipe = globalPipeObjects[leftNeighborId]
                        boundaryHeight = boundaryElement.height
                        rightHeight = rightPipe.heights[1]
                        leftHeight = leftGlobalPipe.heights[-2]
                        forwG = rightHeight - boundaryHeight
                        prevG = boundaryHeight - leftHeight
                        forwVelocity = rightPipe.velocityMesh[stepNumber - 1][1]
                        Vk = leftGlobalPipe.velocityMesh[stepNumber - 1][-2]
                        Vn = leftGlobalPipe.velocityMesh[stepNumber - 1][-1]
                        prevVelocity = interPolateVelocity(Vn=Vn, Vk=Vk, miniSteps=self.miniSteps,
                                                           miniStep=miniStep)
                        forwPressure = rightPipe.pressureMesh[stepNumber - 1][1]
                        Pk = leftGlobalPipe.pressureMesh[stepNumber - 1][-2]
                        Pn = leftGlobalPipe.pressureMesh[stepNumber - 1][-1]
                        prevPressure = interPolatePressure(Pk=Pk, Pn=Pn, miniSteps=self.miniSteps,
                                                           miniStep=miniStep)
                        forwDiameter = rightPipe.innerDiameter
                        prevDiameter = leftGlobalPipe.innerDiameter

                        forwPipeResistance = rightPipe.calcResistance(velocity=forwVelocity)
                        prevPipeResistance = leftGlobalPipe.calcResistance(velocity=prevVelocity)

                        forwKCorrection = rightPipe.kCorrection
                        prevKCorrection = leftGlobalPipe.kCorrection





                    elif rightNeighborId not in list(self.pipeObjects.keys()):
                        rightGlobalPipe = globalPipeObjects[rightNeighborId]

                        boundaryHeight = boundaryElement.height
                        rightHeight = rightGlobalPipe.heights[1]
                        leftHeight = leftPipe.heights[-2]
                        forwG = rightHeight - boundaryHeight
                        prevG = boundaryHeight - leftHeight

                        prevVelocity = leftPipe.velocityMesh[stepNumber - 1][-2]
                        Vk = rightGlobalPipe.velocityMesh[stepNumber - 1][1]
                        Vn = rightGlobalPipe.velocityMesh[stepNumber - 1][0]
                        forwVelocity = interPolateVelocity(Vn=Vn, Vk=Vk, miniSteps=self.miniSteps,
                                                           miniStep=miniStep)
                        prevPressure = leftPipe.pressureMesh[stepNumber - 1][1]
                        Pk = rightGlobalPipe.pressureMesh[stepNumber - 1][1]
                        Pn = rightGlobalPipe.pressureMesh[stepNumber - 1][0]
                        forwPressure = interPolatePressure(Pk=Pk, Pn=Pn, miniSteps=self.miniSteps,
                                                           miniStep=miniStep)
                        forwDiameter = rightGlobalPipe.innerDiameter
                        prevDiameter = leftPipe.innerDiameter
                        forwPipeResistance = rightGlobalPipe.calcResistance(velocity=forwVelocity)
                        prevPipeResistance = leftPipe.calcResistance(velocity=prevVelocity)

                        forwKCorrection = rightGlobalPipe.kCorrection
                        prevKCorrection = leftPipe.kCorrection
                    else:
                        forwVelocity = rightPipe.velocityMesh[stepNumber - 1][1]
                        prevVelocity = leftPipe.velocityMesh[stepNumber - 1][-2]
                        forwPressure = rightPipe.pressureMesh[stepNumber - 1][1]
                        prevPressure = leftPipe.pressureMesh[stepNumber - 1][-2]
                        forwDiameter = rightPipe.innerDiameter
                        prevDiameter = leftPipe.innerDiameter
                        forwPipeResistance = rightPipe.calcResistance(velocity=forwVelocity)
                        prevPipeResistance = leftPipe.calcResistance(velocity=prevVelocity)
                        forwKCorrection = rightPipe.kCorrection
                        prevKCorrection = leftPipe.kCorrection

                        try:
                            boundaryHeight = boundaryElement.height
                            rightHeight = rightPipe.heights[1]
                            leftHeight = leftPipe.heights[-2]
                            forwG = rightHeight - boundaryHeight
                            prevG = boundaryHeight - leftHeight
                        except Exception as e:
                            print(e)
                            forwG = 0
                            prevG = 0

                    boundaryLeftPressure, boundaryRightPressure, boundaryLeftVelocity, boundaryRightVelocity = \
                        boundaryElement.calcBoundaries(tau=self.tau, soundSpeed=soundSpeed, density=density,
                                                       forwG=forwG,
                                                       prevG=prevG,
                                                       forwVelocity=forwVelocity, prevVelocity=prevVelocity,
                                                       forwPressure=forwPressure, prevPressure=prevPressure,
                                                       forwPipeResistance=forwPipeResistance,
                                                       prevPipeResistance=prevPipeResistance,
                                                       forwKCorrection=forwKCorrection,
                                                       prevKCorrection=prevKCorrection, prevDiameter=prevDiameter,
                                                       forwDiameter=forwDiameter)

                    conditions = {
                        'pressure': [boundaryLeftPressure, boundaryRightPressure],
                        'velocity': [boundaryLeftVelocity, boundaryRightVelocity]
                    }

                self.boundaryConditions[boundaryElement.id] = conditions
                if miniStep == self.miniSteps:
                    if leftPipe not in self.pipeObjects.values():
                        boundaryLeftP = boundaryLeftPressure
                        boundaryLeftV = boundaryLeftVelocity
                    if rightPipe not in self.pipeObjects.values():
                        boundaryRightP = boundaryRightPressure
                        boundaryRightV = boundaryRightVelocity

            for pipe in self.pipeObjects.values():
                leftNeighbor = pipe.neighbors[0]
                rightNeighbor = pipe.neighbors[1]
                boundaryLeftVelocity, boundaryLeftPressure, boundaryRightVelocity, boundaryRightPressure = \
                    getBoundaryConditions(leftNeighbor=leftNeighbor, rightNeighbor=rightNeighbor,
                                          pipeId=pipe.id, boundaryObjects=self.boundaryObjects,
                                          boundaryDict=self.boundaryConditions)
                pipe.calcStep(soundSpeed=soundSpeed, density=density, boundaryLeftVelocity=boundaryLeftVelocity,
                              boundaryRightVelocity=boundaryRightVelocity, boundaryLeftPressure=boundaryLeftPressure,
                              boundaryRightPressure=boundaryRightPressure, stepNumber=stepNumber, tau=self.tau)
                if miniStep != self.miniSteps:
                    pipe.pressureMesh[stepNumber - 1] = pipe.pressureMesh[stepNumber]
                    pipe.velocityMesh[stepNumber - 1] = pipe.velocityMesh[stepNumber]
        return boundaryLeftP, boundaryRightP, boundaryLeftV, boundaryRightV
