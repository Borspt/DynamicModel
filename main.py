import json
from Elements import *

JSON = "schemetest (3) (pump).json"

with open(JSON, "r", encoding="utf-8") as load_file:
    test_json = json.load(load_file)

elements_dict = {}

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


def init_objects():
    all_objects = []
    pipe_objects = []
    boundary_objects = []
    for blockValve in blockValves:
        elements_dict[blockValve['id']] = 'blockValve'
        object = BlockValve(elementid=blockValve['id'], status=blockValve['status']['value'])
        all_objects.append(object)
        boundary_objects.append(object)

    for branch in branches:
        elements_dict[branch['id']] = 'branch'
        object = Branch(elementId=branch['id'])
        all_objects.append(object)
        boundary_objects.append(object)


    for tank in tanks:
        tank['status'] = {'value': 'open'}
        elements_dict[tank['id']] = 'tank'
        density = tank.get('density', None)

        if tank['position'][0]['neighborsId'][0] == -1:
            orientation = True
        else:
            orientation = False

        object = Tank(elementId=tank['id'], height=tank['height'], inHeight=tank['inHeight'],
                      start=orientation, density=density)
        all_objects.append(object)
        boundary_objects.append(object)


    for fps in flowPressureSetters:
        fps['status'] = {'value': 'open'}
        elements_dict[fps['id']] = 'fps'
        object = FlowPressureSetter(elementid=fps['id'], height=fps['height'], density=fps['density'],
                                    inheight=fps['inHeight'])
        all_objects.append(object)
        boundary_objects.append(object)


    for fs in FilterStrainers:
        elements_dict[fs['id']] = 'fs'
        object = FilterStrainer(elementid=fs['id'])
        all_objects.append(object)
        boundary_objects.append(object)


    for pipe in pipes:
        elements_dict[pipe['id']] = 'pipe'
        if len(pipe['profile']) > 0:
            start_point = pipe['profile'][0].get('distance', 0)
            end_point = pipe['profile'][-1].get('distance', 0)
            start_height = pipe['profile'][0].get('height', 0)
            end_height = pipe['profile'][-1].get('height', 0)
            profile = pipe['profile']
            # start_height = 100
            # end_height = 100

            innerDiameter = pipe['profile'][0].get('innerDiameter', None)
            length = end_point - start_point
            assert length > 0, 'Pipe length must be more than 0'
            assert innerDiameter is not None, 'innerDiameter is None'
        else:
            length = 0
            start_height = 0
            end_height = 0
            innerDiameter = 1.0

        object = Pipe(elementid=pipe['id'], innerDiameter=innerDiameter, length=length, profile=profile,
                      start_height=start_height, end_height=end_height)
        all_objects.append(object)
        pipe_objects.append(object)


    rotor_objects = []
    for rotor in rotors:
        elements_dict[rotor['id']] = 'rotor'
        object = Rotor(elementId=rotor['id'], nominal_frequency=rotor['rotorFrequency']['nominal'],
                       real_frequency=rotor['rotorFrequency']['real'], user_choice=rotor['user_choice'])
        all_objects.append(object)
        rotor_objects.append(object)

    for pump in pumps:
        elements_dict[pump['id']] = 'pump'
        object = Pump(elementId=pump['id'], rotorId=pump['rotorId'])
        object.rotor_objects = rotor_objects
        all_objects.append(object)
        boundary_objects.append(object)

    # for checkValve in checkValves:
    #     elements_dict[checkValve['id']] = 'checkValve'
    #     object = CheckValve(elementId=checkValve['id'])
    #     all_objects.append(object)
    #     boundary_objects.append(object)






    return all_objects, pipe_objects, boundary_objects


all_objects, pipe_objects, boundary_objects = init_objects()


def find_object(target_id, objects):
    for _object in objects:
        if _object.id == target_id:
            return _object
    return None



pipe = find_object(129, pipe_objects)
pipe.calcMesh(soundSpeed=1000, density=850, boundaryLeft= 2*g*10000, boundaryRight= 2*g*10000, tau= 0.1, processTime=1000, showgraph=True)







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
