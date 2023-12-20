import json
from Elements import *
from matplotlib import pyplot as plt




def init_objects():
    all_objects = []
    pipe_objects = {}
    boundaryObjects = {}
    blockValves = test_json['topology'].get('blockValves', [])
    branches = test_json['topology'].get('branches', [])
    pipes = test_json['topology'].get('pipes', [])
    pumps = test_json['topology'].get('pumps', [])
    tanks = test_json['topology'].get('tanks', [])
    rotors = test_json['topology'].get('rotors', [])
    flowPressureSetters = test_json['topology']['flowPressureSetters']
    # plugs = test_json['topology']['plugs']
    FilterStrainers = test_json['topology']['FilterStrainers']
    # checkValves = test_json['topology']['checkValves']
    Neighbor = 'neighborsId'
    for blockValve in blockValves:
        elementId = blockValve['id']
        element = BlockValve(elementId=elementId, status=blockValve['status']['value'],
                             neighbors=blockValve['position'][0].get(
                                 Neighbor, None))
        all_objects.append(element)
        boundaryObjects[elementId] = element

    for branch in branches:
        elementId = branch['id']
        element = Branch(elementId=branch['id'], neighbors=branch['position'][0].get(Neighbor, None))
        all_objects.append(element)
        boundaryObjects[elementId] = element

    for tank in tanks:
        tank['status'] = {'value': 'open'}
        elementId = tank['id']
        density = tank.get('density', None)

        if tank['position'][0][(Neighbor)][0] == -1:
            orientation = True
        else:
            orientation = False

        element = Tank(elementId=tank['id'], height=tank['height'], inHeight=tank['inHeight'],
                       start=orientation, density=density, neighbors=tank['position'][0][(Neighbor)])
        all_objects.append(element)
        boundaryObjects[elementId] = element

    for fps in flowPressureSetters:
        fps['status'] = {'value': 'open'}
        elementId = fps['id']
        if fps['position'][0][(Neighbor)][0] == -1:
            orientation = True
        else:
            orientation = False

        element = FlowPressureSetter(elementId=fps['id'], height=fps['height'], start=orientation,
                                     density=fps['density'],
                                     inHeight=fps['inHeight'], neighbors=fps['position'][0].get(Neighbor, None))
        all_objects.append(element)
        boundaryObjects[elementId] = element

    for fs in FilterStrainers:
        elementId = fs['id']
        element = FilterStrainer(elementiId=fs['id'], neighbors=fs['position'][0].get(Neighbor, None))
        all_objects.append(element)
        boundaryObjects[elementId] = element

    rotor_objects = {}
    for rotor in rotors:
        elementId = rotor['id']
        element = Rotor(elementId=rotor['id'], nominal_frequency=rotor['rotorFrequency']['nominal'],
                        real_frequency=rotor['rotorFrequency']['real'], user_choice=rotor['user_choice'])
        all_objects.append(element)

    for pump in pumps:
        elementId = pump['id']
        if pump['useConstHead'] == 1:
            constHead = True
            deltaHead = pump['deltaHead']
        else:
            constHead = False
            deltaHead = None
        height = pump['height']
        element = Pump(elementId=pump['id'], rotorId=pump['rotorId'], useConstHead=constHead, deltaHead=deltaHead,
                       height=height, neighbors=pump['position'][0].get(Neighbor, None))
        element.rotor_objects = rotor_objects
        all_objects.append(element)
        boundaryObjects[elementId] = element

    for pipe in pipes:
        elementId = pipe['id']
        neighbors = pipe['position'][0].get(Neighbor, None)
        if len(pipe['profile']) > 0:
            start_point = pipe['profile'][0].get('distance', 0)
            end_point = pipe['profile'][-1].get('distance', 0)
            start_height = pipe['profile'][0].get('height', 0)
            end_height = pipe['profile'][-1].get('height', 0)
            profile = pipe['profile']

            innerDiameter = pipe['profile'][0].get('innerDiameter', None)
            length = end_point - start_point
            assert length > 0, 'Pipe length must be more than 0'
            assert innerDiameter is not None, 'innerDiameter is None'
        else:
            length = 0
            start_height = boundaryObjects[neighbors[0]].height
            end_height = boundaryObjects[neighbors[1]].height
            profile = None
            innerDiameter = 1.0

        element = Pipe(elementId=pipe['id'], innerDiameter=innerDiameter, length=length, profile=profile,
                       start_height=start_height, end_height=end_height, neighbors=neighbors)
        all_objects.append(element)
        pipe_objects[elementId] = element

    # for checkValve in checkValves:
    #     elements_dict[checkValve['id']] = 'checkValve'
    #     element = CheckValve(elementId=checkValve['id'])
    #     all_objects.append(element)
    #     boundary_objects.append(element)

    return all_objects, pipe_objects, boundaryObjects





def calcInitValues(boundaryObjects, initialVelocity, initialPressure):
    boundaryConditions = {}
    for element in boundaryObjects.values():
        if isinstance(element, (FlowPressureSetter, Tank)):
            if element.start:
                elementConditions = {
                    'velocity': [None, initialVelocity],
                    'pressure': [None, element.deltaP(density=element.density)]
                }
            else:
                try:
                    elementConditions = {
                        'velocity': [initialVelocity, None],

                        'pressure': [element.deltaP(density=element.density), None]}
                except:
                    print(element.__dict__)


        else:
            elementConditions = {
                'velocity': [initialVelocity, initialVelocity],
                'pressure': [initialPressure, initialPressure]
            }
        print(element)
        boundaryConditions[element.id] = elementConditions

    return boundaryConditions


def find_object(target_id, objects):
    for _object in objects.values:
        if _object.id == target_id:
            return _object
    return None


def getBoundaryConditions(leftNeighbor, rightNeighbor, pipeId):
    boundaryLeft = boundaryConditions.get(leftNeighbor, None)
    boundaryRight = boundaryConditions.get(rightNeighbor, None)

    assert boundaryLeft is not None, f'Pipe id{pipeId} left neighbor not in boundaryConditions'
    assert boundaryRight is not None, f'Pipe id{pipeId} right neighbor not in boundaryConditions'

    boundaryLeftVelocity = boundaryLeft.get('velocity', None)[1]
    boundaryLeftPressure = boundaryLeft.get('pressure', None)[1]

    assert boundaryLeftVelocity is not None, f'Pipe id{pipeId} left velocity is None'
    assert boundaryLeftPressure is not None, f'Pipe id{pipeId} left pressure is None'

    boundaryRightVelocity = boundaryRight.get('velocity', None)[0]
    boundaryRightPressure = boundaryRight.get('pressure', None)[0]

    assert boundaryRightVelocity is not None, f'Pipe id{pipeId} right velocity is None'
    assert boundaryRightPressure is not None, f'Pipe id{pipeId} right pressure is None'

    return boundaryLeftVelocity, boundaryLeftPressure, boundaryRightVelocity, boundaryRightPressure


JSON = "schemes/schemetest (3) (2tanks).json"

with open(JSON, "r", encoding="utf-8") as load_file:
    test_json = json.load(load_file)

elements_dict = {}

PROCESSTIME = 150
TAU = 0.1
SOUNDSPEED = 1000
DENSITY = 850
INITIAL_VELOCITY = 0.01
INITIAL_PRESSURE = g * 10000
SHOW_GRAPH = True
USER_CHOICE = [123]

all_objects, pipeObjects, boundaryObjects = init_objects()

steps = math.ceil(PROCESSTIME / TAU)

boundaryConditions = calcInitValues(boundaryObjects=boundaryObjects, initialVelocity=INITIAL_VELOCITY,
                                    initialPressure=INITIAL_PRESSURE)

for pipe in pipeObjects.values():
    pipe.initMesh(soundSpeed=SOUNDSPEED, steps=steps, tau=1)

    pipe.velocityMesh[0].fill(0.01)
    pipe.pressureMesh[0].fill(g * 10000)

    boundaryLeftVelocity, boundaryLeftPressure, boundaryRightVelocity, boundaryRightPressure = \
        getBoundaryConditions(leftNeighbor=pipe.neighbors[0], rightNeighbor=pipe.neighbors[1], pipeId=pipe.id)

    pipe.pressureMesh[0][0], pipe.pressureMesh[0][-1] = boundaryLeftPressure, boundaryRightPressure
    pipe.velocityMesh[0][0], pipe.velocityMesh[0][-1] = boundaryLeftVelocity, boundaryRightVelocity

    '''
    Заполнение
    граничных
    условий
    '''

for stepNumber in range(1, steps + 1):
    for boundaryElement in boundaryObjects.values():
        if isinstance(boundaryElement, (FlowPressureSetter, Tank)):
            if boundaryElement.start:
                neighborId = boundaryElement.neighbors[1]
                neighborPoint = 1
            else:
                neighborId = boundaryElement.neighbors[0]
                neighborPoint = -2
            pipe = pipeObjects[neighborId]
            boundaryLeftPressure, boundaryRightPressure,\
                boundaryLeftVelocity, boundaryRightVelocity = \
                boundaryElement.calcBoundaries(tau=TAU, soundSpeed=SOUNDSPEED, pipeHeight=pipe.heights[neighborPoint],
                                               pipeVelocity=pipe.velocityMesh[stepNumber - 1][neighborPoint],
                                               pipePressure=pipe.pressureMesh[stepNumber - 1][neighborPoint],
                                               pipeResistance=pipe.calcResistance(
                                                   velocity=(pipe.velocityMesh[stepNumber - 1][neighborPoint])),
                                               pipeKCorrection=pipe.kCorrection)
        else:
            leftNeighborId = boundaryElement.neighbors[0]
            rightNeighborId = boundaryElement.neighbors[1]
            leftPipe = pipeObjects[leftNeighborId]
            rightPipe = pipeObjects[rightNeighborId]
            forwVelocity = rightPipe.velocityMesh[stepNumber - 1][0]
            prevVelocity = leftPipe.velocityMesh[stepNumber - 1][-1]
            forwPressure = rightPipe.pressureMesh[stepNumber - 1][0]
            prevPressure = leftPipe.pressureMesh[stepNumber - 1][-1]
            forwDiameter = leftPipe.innerDiameter
            prevDiameter = leftPipe.innerDiameter
            boundaryHeight = boundaryElement.height
            rightHeight = rightPipe.start_height
            leftHeight = leftPipe.end_height

            boundaryLeftPressure, boundaryRightPressure, boundaryLeftVelocity, boundaryRightVelocity = \
                boundaryElement.calcBoundaries(tau=TAU, soundSpeed=SOUNDSPEED, density=DENSITY,
                                               forwG=(boundaryHeight - rightHeight)/2,
                                               prevG=(leftHeight - boundaryHeight)/2,
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
            leftNeighbor = pipe.neighbors[0]
            rightNeighbor = pipe.neighbors[1]
            boundaryLeftVelocity, boundaryLeftPressure, boundaryRightVelocity, boundaryRightPressure = \
                getBoundaryConditions(leftNeighbor=leftNeighbor, rightNeighbor=rightNeighbor, pipeId=pipe.id)
            pipe.calcStep(soundSpeed=SOUNDSPEED, density=DENSITY, boundaryLeftVelocity=boundaryLeftVelocity,
                          boundaryRightVelocity=boundaryRightVelocity, boundaryLeftPressure=boundaryLeftPressure,
                          boundaryRightPressure=boundaryRightPressure, stepNumber=stepNumber, tau=TAU)

if SHOW_GRAPH:
    plt.ion()
    fig = plt.figure()
    time = 0
    x = []
    y = []
    x_end = []
    y_end = []
    for elementId in USER_CHOICE:
        pipe = pipeObjects[elementId]
        for i in range(len(pipe.profile)):  ## Почему тут берутся данные по профилю, а не сетке?
            # - строится красный профиль высот и отмечаются концы труб

            el = pipe.profile[i]
            if i == 0 or (i == len(pipe.profile) - 1):
                x_end.append(el['distance'])
                y_end.append(el['height'])
            x.append(el['distance'])
            y.append(el['height'])
    print(f'Y = {y}')

    for step in range(1, stepNumber + 1):
        time += TAU
        dotX = np.array([])
        ph = np.array([])
        for elementId in USER_CHOICE:
            pipe = pipeObjects[elementId]
            additionalX = np.array(pipe.meshX) / 1000 + pipe.profile[0]['distance']
            dotX = np.append(dotX, additionalX)

            for i in range(len(pipe.pressureMesh[0])):
                # ph.append((self.pressureMesh[step-1][i] + density*g*heights[i])/density/g)
                additionalPh = (pipe.pressureMesh[step - 1][i] + DENSITY * g * pipe.heights[i]) / DENSITY / g
                ph = np.append(ph, additionalPh)
        # print(ph)

        plt.clf()
        plt.plot(dotX, ph, 'blue', marker='o', label='Гидроуклон')
        plt.plot(x_end, y_end, 'black', marker='D', label='Концы')
        plt.plot(x, y, color='red', ls='--', label='Высота пролегания')
        # plt.scatter(dotX, self.heights, marker='o')
        step_text = f'Шаг: {step} из {stepNumber + 1}'
        plt.text(0.95, 0.95, step_text, transform=plt.gca().transAxes, fontsize=10,
                 verticalalignment='top', horizontalalignment='right')
        plt.title(f'Обновляемый график {time:.1f} с')
        plt.xlabel('Дистанция, км', fontsize=30)
        plt.ylabel('Напор, м', fontsize=30)
        plt.yticks(fontsize=20)
        plt.xticks(fontsize=20)
        plt.pause(TAU)
        # plt.show()
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
