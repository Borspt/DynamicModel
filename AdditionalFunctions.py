from Elements import *


def init_objects(test_json, tau):
    all_objects = []
    pipeObjects = {}
    technoObjects = {}
    boundaryObjects = {}
    branchIdList = []

    blockValves = test_json['topology'].get('blockValves', [])
    branches = test_json['topology'].get('branches', [])
    pipes = test_json['topology'].get('pipes', [])
    pumps = test_json['topology'].get('pumps', [])
    tanks = test_json['topology'].get('tanks', [])
    rotors = test_json['topology'].get('rotors', [])
    flowPressureSetters = test_json['topology']['flowPressureSetters']
    # plugs = test_json['topology']['plugs']
    FilterStrainers = test_json['topology']['FilterStrainers']
    checkValves = test_json['topology']['checkValves']
    if 'assembly' in test_json.keys():
        technologicalObjects = test_json['assembly']['technologicalObjects']
    else:
        technologicalObjects = None
        boundaryTechnoIds = []
        pipeTechnoIds = []
    Neighbor = 'neighborsId'

    for blockValve in blockValves:
        elementId = blockValve['id']
        height = blockValve['height']
        element = BlockValve(elementId=elementId, height=height, status=blockValve['status']['value'],
                             neighbors=blockValve['position'][0].get(
                                 Neighbor, None))
        all_objects.append(element)
        boundaryObjects[elementId] = element

    for checkValve in checkValves:
        elementId = checkValve['id']
        height = checkValve['height']
        element = CheckValve(elementId=elementId, height=height,
                             neighbors=checkValve['position'][0].get(Neighbor, None))
        all_objects.append(element)
        boundaryObjects[elementId] = element

    for branch in branches:
        elementId = branch['id']
        height = branch['height']
        neighborsIn = branch['position'][0].get('neighborsIn', None)
        neighborsOut = branch['position'][0].get('neighborsOut', None)
        element = Branch(elementId=branch['id'], height=height, neighborsIn=neighborsIn, neighborsOut=neighborsOut)
        all_objects.append(element)
        boundaryObjects[elementId] = element
        branchIdList.append(elementId)

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

        profileExist = pipe.get('profile', None)
        if profileExist is None:
            lenProfile = 0
        else:
            lenProfile = len(pipe['profile'])
        if lenProfile > 0:
            pipeTechnological = False
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
            pipeTechnological = True
            length = 0
            start_height = boundaryObjects[neighbors[0]].height
            end_height = boundaryObjects[neighbors[1]].height
            profile = None
            innerDiameter = 1.0
        pipeId = pipe['id']
        if len(pipe['position']) > 2:
            technoNeighbors = pipe['position'][2].get(Neighbor, None)
        else:
            technoNeighbors = None
        element = Pipe(elementId=pipeId, innerDiameter=innerDiameter, length=length, profile=profile,
                       start_height=start_height, end_height=end_height, neighbors=neighbors,
                       technoNeighbors=technoNeighbors, isTechnological=pipeTechnological)
        all_objects.append(element)
        pipeObjects[elementId] = element
    if technologicalObjects is not None:
        for technoObject in technologicalObjects:
            technoId = technoObject['id']
            height = technoObject['height']
            boundaryTechnoObjects = {}
            pipeTechnoObjects = {}
            boundaryTechnoIds = technoObject['blockValvesId'] + technoObject['pumpsId'] \
                                + technoObject['checkValvesId'] + technoObject['tanksId'] + technoObject['branchesId']
            for elementId in boundaryTechnoIds:
                boundaryTechnoObjects[elementId] = boundaryObjects[elementId]
            pipeTechnoIds = technoObject['pipesId']
            for elementId in pipeTechnoIds:
                pipeTechnoObjects[elementId] = pipeObjects[elementId]

        # miniSteps = len(pipeTechnoIds)
            miniSteps = 100
            element = TechnologicalObject(elementId=technoId, boundaryObjects=boundaryTechnoObjects,
                                          pipeObjects=pipeTechnoObjects, miniSteps=miniSteps, tau=tau / miniSteps,
                                          height=height)
            technoObjects[technoId] = element
            boundaryObjects[technoId] = element

    return all_objects, pipeObjects, boundaryObjects, branchIdList, technoObjects, boundaryTechnoIds, pipeTechnoIds


def calcInitValues(boundaryObjects, initialVelocity, initialPressure, density):
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

                        'pressure': [-element.deltaP(density=element.density), None]}
                except Exception as e:
                    print(e)
                    print(element.__dict__)

        elif isinstance(element, Branch):
            heightAdjustedPressure = initialPressure - element.height * density * g
            elementConditions = {
                'velocity': [initialVelocity, initialVelocity, initialVelocity],
                'pressure': [heightAdjustedPressure, heightAdjustedPressure, heightAdjustedPressure]
            }
        else:
            heightAdjustedPressure = initialPressure - element.height * density * g
            elementConditions = {
                'velocity': [initialVelocity, initialVelocity],
                'pressure': [heightAdjustedPressure, heightAdjustedPressure]
            }
        # print(element)
        boundaryConditions[element.id] = elementConditions

    return boundaryConditions


def find_object(target_id, objects):
    for _object in objects.values:
        if _object.id == target_id:
            return _object
    return None



