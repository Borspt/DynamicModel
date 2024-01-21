import random
import json


def make_profile(start_point, distance, start_height, end_height, sigma, points=10):
    profile = []
    step_length = distance / points
    step_height = (end_height - start_height) / points

    for num in range(points + 1):
        if num == 0 or num == points:
            profile.append(
                {"distance": start_point + num * step_length,
                 "height": start_height + num * step_height,
                 "innerDiameter": 1,

                 "roughness": 0.0002,
                 "maxPressure": 2000000}
            )
        else:
            height = start_height + num * step_height
            random_height = random.normalvariate(mu=height, sigma=sigma)
            profile.append(
                {"distance": start_point + num * step_length,
                 "height": random_height,
                 "innerDiameter": 1,

                 "roughness": 0.0002,
                 "maxPressure": 2000000}
            )
    return profile

profile = make_profile(start_point=105, distance=20, start_height=150, end_height=150, sigma=5)
print(profile)
with open('profiles/profile.json', "w", encoding="utf-8") as post_file:
    json.dump(profile, post_file, indent=4, ensure_ascii=False)