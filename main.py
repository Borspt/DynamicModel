import json
from Elements import *

JSON = "schemes/schemetest (3) (pump).json"

with open(JSON, "r", encoding="utf-8") as load_file:
    test_json = json.load(load_file)

elements_dict = {}




def init_objects():
    all_objects = []
    pipe_objects = {}
    boundary_objects = {}
    blockValves = test_json['topology']['blockValves']
    branches = test_json['topology']['branches']
    pipes = test_json['topology']['pipes']
    pumps = test_json['topology']['pumps']
    tanks = test_json['topology']['tanks']
    rotors = test_json['topology']['rotors']
    flowPressureSetters = test_json['topology']['flowPressureSetters']
    # plugs = test_json['topology']['plugs']
    FilterStrainers = test_json['topology']['FilterStrainers']
    # checkValves = test_json['topology']['checkValves']
    for blockValve in blockValves:
        elementId = blockValve['id']
        element = BlockValve(elementId=elementId, status=blockValve['status']['value'], neighbors=blockValve['profile'][0].get('neighbors', None))
        all_objects.append(element)
        boundary_objects[elementId] = element

    for branch in branches:
        elementId = branch['id']
        element = Branch(elementId=branch['id'], neighbors=branch['profile'][0].get('neighbors', None))
        all_objects.append(element)
        boundary_objects[elementId] = element


    for tank in tanks:
        tank['status'] = {'value': 'open'}
        elementId = tank['id']
        density = tank.get('density', None)

        if tank['position'][0]['neighborsId'][0] == -1:
            orientation = True
        else:
            orientation = False

        element = Tank(elementId=tank['id'], height=tank['height'], inHeight=tank['inHeight'],
                      start=orientation, density=density, neighbors=tank['profile'][0].get('neighbors', None))
        all_objects.append(element)
        boundary_objects[elementId] = element


    for fps in flowPressureSetters:
        fps['status'] = {'value': 'open'}
        elementId = fps['id']
        if fps['position'][0]['neighborsId'][0] == -1:
            orientation = True
        else:
            orientation = False

        element = FlowPressureSetter(elementId=fps['id'], height=fps['height'], start=orientation, density=fps['density'],
                                    inHeight=fps['inHeight'], neighbors=fps['profile'][0].get('neighbors', None))
        all_objects.append(element)
        boundary_objects[elementId] = element


    for fs in FilterStrainers:
        elementId = fs['id']
        element = FilterStrainer(elementiId=fs['id'], neighbors=fs['profile'][0].get('neighbors', None))
        all_objects.append(element)
        boundary_objects[elementId] = element





    rotor_objects = {}
    for rotor in rotors:
        elementId = rotor['id']
        element = Rotor(elementId=rotor['id'], nominal_frequency=rotor['rotorFrequency']['nominal'],
                       real_frequency=rotor['rotorFrequency']['real'], user_choice=rotor['user_choice'])
        all_objects.append(element)
        boundary_objects[elementId] = element

    for pump in pumps:
        elementId = pump['id']
        element = Pump(elementId=pump['id'], rotorId=pump['rotorId'], neighbors=pump['profile'][0].get('neighbors', None))
        element.rotor_objects = rotor_objects
        all_objects.append(element)
        boundary_objects[elementId] = element

    for pipe in pipes:
        elementId = pipe['id']
        neighbors = pipe['profile'][0].get('neighbors', None)
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
            start_height = boundary_objects[neighbors[0]].height
            end_height = boundary_objects[neighbors[1]].height
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






    return all_objects, pipe_objects, boundary_objects


all_objects, pipeObjects, boundary_objects = init_objects()


def calcInitValues(boundary_objects, INITIAL_VELOCITY, INITIAL_PRESSURE):
    boundaryConditions = {}
    for object in boundary_objects:
        if isinstance(object, (FlowPressureSetter, Tank)):
            if object.start:
                elementConditions = {
                    'velocity' : [None, INITIAL_VELOCITY],
                    'pressure' : [None,object.deltaP(density=object.density)]
                }
            else:
                elementConditions = {
                    'velocity': [INITIAL_VELOCITY, None],
                    'pressure': [object.deltaP(density=object.density), None]
                }
        else:
            elementConditions = {
                'velocity': [INITIAL_VELOCITY, INITIAL_VELOCITY],
                'pressure': [INITIAL_PRESSURE, INITIAL_PRESSURE]
            }
        boundaryConditions[object.id] = elementConditions

    return boundaryConditions



def find_object(target_id, objects):
    for _object in objects:
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





PROCESSTIME = 15
TAU = 0.1
SOUNDSPEED = 1000
DENSITY = 850
INITIAL_VELOCITY = 0.01
INITIAL_PRESSURE = g * 10000


steps = math.ceil(PROCESSTIME / TAU)

boundaryConditions = calcInitValues(boundary_objects=boundary_objects, initialPressure=INITIAL_PRESSURE,
                                    initialVelocity=INITIAL_VELOCITY)



for pipe in pipeObjects:

    pipe.initMesh(soundSpeed = SOUNDSPEED,steps=steps, tau=1)

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

    #####



    # for i in range(self.pressureMesh.shape[0]):
    #     self.pressureMesh[i][0] = boundaryLeft
    #     self.pressureMesh[i][-1] = boundaryRight






for stepNumber in range(1, steps + 1):
    for boundaryElement in boundary_objects:
        if isinstance(boundaryElement, (FlowPressureSetter, Tank)):
            if boundaryElement.start:
                neighborId = boundaryElement.neighbors[1]
                neighborPoint = 1
            else:
                neighborId = boundaryElement.neighbors[0]
                neighborPoint = -2
            pipe = find_object(neighborId, pipeObjects)
            pipe = pipeObjects[neighborId]
            boundaryLeftPressure, boundaryRightPressure, boundaryLeftVelocity, boundaryRightVelocity =\
                boundaryElement.calcBoundaries(tau=TAU, soundSpeed=SOUNDSPEED, pipeHeight = pipe.heights[neighborPoint],
                                               pipeVelocity = pipe.velocityMesh[stepNumber-1][neighborPoint],
                                               pipePressure = pipe.pressureMesh[stepNumber-1][neighborPoint],
                                               pipeResistance = pipe.calcResistance(
                                                   velocity=(pipe.velocityMesh[stepNumber - 1][neighborPoint])),
                                                pipeKCorrection = pipe.kCorrection)
        else:
            leftNeighborId = boundaryElement.neighbors[0]
            rightNeighborId = boundaryElement.neighbors[1]
            leftPipe = pipeObjects[leftNeighborId]
            rightPipe = pipeObjects[rightNeighborId]
            forwVelocity = rightPipe.velocityMesh[stepNumber - 1][0]
            prevVelocity = leftPipe.velocityMesh[stepNumber - 1][-1]
            forwPressure = rightPipe.velocityMesh[stepNumber - 1][0]
            prevPressure = leftPipe.velocityMesh[stepNumber - 1][-1]
            boundaryLeftPressure, boundaryRightPressure, boundaryLeftVelocity, boundaryRightVelocity = \
                boundaryElement.calcBoundaries(tau=TAU, soundSpeed=SOUNDSPEED, density = DENSITY,
                                               forwG = rightPipe.start_height, prevG = leftPipe.end_height,
                                               forwVelocity = forwVelocity, prevVelocity = prevVelocity,
                                               forwPressure = forwPressure, prevPressure = prevPressure,
                                               forwResistance = rightPipe.calcResistance(velocity = forwVelocity),
                                               prevResistance = leftPipe.calcResistance(velocity = prevVelocity),
                                               forwKCorrection = rightPipe.kCorrection,
                                               prevKCorrection = leftPipe.kCorrection)

        conditions = {
            'pressure' : [boundaryLeftPressure, boundaryRightPressure],
            'velocity' : [boundaryLeftVelocity, boundaryRightVelocity]
            }

        boundaryConditions[boundaryElement.id] = conditions

        for pipe in pipeObjects:
            leftNeighbor = pipe.neighbors[0]
            rightNeighbor = pipe.neighbors[1]
            boundaryLeftVelocity, boundaryLeftPressure, boundaryRightVelocity, boundaryRightPressure = \
                getBoundaryConditions(leftNeighbor=leftNeighbor, rightNeighbor=rightNeighbor, pipeId=pipe.id)
            pipe.calcStep(soundSpeed = SOUNDSPEED, density = DENSITY, boundaryLeftVelocity=boundaryLeftVelocity,
                          boundaryRightVelocity=boundaryRightVelocity, boundaryLeftPressure=boundaryLeftPressure,
                          boundaryRightPressure=boundaryRightPressure,stepNumber = stepNumber, tau=TAU, showgraph=False)










print(pipe.__dict__)
print(pipe.profile)
print(pipe.heights)
print(pipe.meshX)
print(len(pipe.heights))
print(len(pipe.meshX))

print('_________________________________________________________________')
print('Velocity Mesh:')
print(pipe.velocityMesh)
print('_________________________________________________________________')
print(pipe.pressureMesh)
