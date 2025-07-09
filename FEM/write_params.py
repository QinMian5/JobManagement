# Author: Mian Qin
# Date Created: 2025/5/20
import json


parameters = {
    "flat": {
        "45": [70, 40],
        "60": [65, 50],
        "70": [55, 50],
        "80": [55, 55],
        "90": [55, 60],
        "120": [45, 70],
        "150": [45, 80],
    },
    "pillar": {
        "45": [70 + 10, 40],
        "60": [65 + 10, 50],
        "70": [55 + 10, 50],
        "80": [55 + 10, 55],
        "90": [55 + 10, 60],
        "120": [45 + 10, 70],
        "150": [45 + 10, 80],
    }
}


def main():
    # system, theta, job_name
    jobs = [
        ["flat", "60", "r30_to_l3000"],
        ["flat", "90", "r30_to_l3000"],
        ["flat", "120", "r30_to_l3000"],
    ]
    job_params = []
    for job in jobs:
        system, theta, job_name = job
        r_box, z_box = parameters[system][theta]
        job_params.append({
            "SYSTEM": system,
            "THETA": theta,
            "JOB_NAME": job_name,
        })
    with open("job_params.json", 'w') as file:
        json.dump(job_params, file, indent='\t')


if __name__ == "__main__":
    main()
