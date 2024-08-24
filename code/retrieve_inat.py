# RETRIEVE INAT DATA #
# Zane Libke #
# 18 April 2024 #
# Written with copious help from ChatGPT4 #



import pandas as pd
import requests

def get_inaturalist_observation_data(observation_url):
    # Extract the observation ID from the URL
    # Assumes the URL format is something like "https://www.inaturalist.org/observations/123456"
    observation_id = observation_url.split('/')[-1]
    
    # Construct the API request URL
    api_url = f'https://api.inaturalist.org/v1/observations/{observation_id}'
    
    # Send the request to the iNaturalist API
    response = requests.get(api_url)
    
    # Check if the request was successful
    if response.status_code == 200:
        # Parse the JSON response into a dictionary
        data = response.json()
        # Return the data
        return data
    else:
        # Handle errors (e.g., observation not found, API issues)
        return {'error': 'Failed to retrieve data', 'status_code': response.status_code}

def parse_observation_details(data):
    # Initialize a result dictionary to store the details
    result = {}
    try:
        obs = data['results'][0]
        # Extract species name
        result['species_name'] = obs['taxon']['name'] if obs['taxon'] else 'Unknown'
        # Extract coordinates
        if 'location' in obs and obs['location']:
            lat, lng = obs['location'].split(',')
            result['latitude'] = lat.strip()
            result['longitude'] = lng.strip()
        else:
            result['latitude'] = 'NA'
            result['longitude'] = 'NA'
        result['gps_accuracy'] = obs['positional_accuracy'] if 'positional_accuracy' in obs else 'Unknown'
        # Extract time & date observed
        result['date_observed'] = obs['observed_on_details']['date'] if 'observed_on_details' in obs else 'Unknown'
        result['time_observed'] = obs['results']['time_observed_at'] if 'results' in obs else 'Unknown'
        # Extract notes
        result['notes'] = obs['description'] if 'description' in obs else 'No notes'
    except KeyError:
        # Handle missing data in the response
        result = {'error': 'Necessary data missing in the response'}
    
    return result

def process_inaturalist_links(csv_file_path):
    # Read the CSV file
    df = pd.read_csv(csv_file_path)
    
    # Initialize columns for the data to be added
    df['species_name'] = None
    df['latitude'] = None
    df['longitude'] = None
    df['gps_accuracy'] = None
    df['date_observed'] = None
    df['time_observed'] = None
    df['notes'] = None
    
    # Process each link in the DataFrame
    for index, row in df.iterrows():
        if pd.notna(row['inat_link']):
            details = get_inaturalist_observation_data(row['inat_link'])
            details = parse_observation_details(details)
            df.at[index, 'species_name'] = details.get('species_name', 'Error')
            df.at[index, 'longitude'] = details.get('longitude', 'Error')
            df.at[index, 'latitude'] = details.get('latitude', 'Error')
            df.at[index, 'gps_accuracy'] = details.get('gps_accuracy', 'Error')
            df.at[index, 'date_observed'] = details.get('date_observed', 'Error')
            df.at[index, 'time_observed'] = details.get('time_observed', 'Error')
            df.at[index, 'notes'] = details.get('notes', 'Error')
    
    # Save the updated DataFrame to a new CSV file
    df.to_csv('updated_data.csv', index=False)
    return 'updated_data.csv'


#function to retrieve all inat data from a df, create new columns for retrieved data, and return the new updated df
def retrieve_inat_all(samples):
    # Initialize columns for inat data to be added
    samples['species_inat'] = None
    samples['latitude_inat'] = None
    samples['longitude_inat'] = None
    samples['gps-accuracy_inat'] = None
    samples['date-observed_inat'] = None
    samples['time-observed_inat'] = None
    samples['notes_inat'] = None

    for idx, row in samples.iterrows():
        #if inat id exists, fetch location data etc.
        if pd.notna(row['inat_link']):
            try:
                details = get_inaturalist_observation_data(row['inat_link'])
                details = parse_observation_details(details)
                samples.at[idx, 'species_inat'] = details.get('species_name', 'Error')
                samples.at[idx, 'longitude_inat'] = details.get('longitude', 'Error')
                samples.at[idx, 'latitude_inat'] = details.get('latitude', 'Error')
                samples.at[idx, 'gps-accuracy_inat'] = details.get('gps_accuracy', 'Error')
                samples.at[idx, 'date-observed_inat'] = details.get('date_observed', 'Error')
                samples.at[idx, 'time-observed_inat'] = details.get('time_observed', 'Error')
                samples.at[idx, 'notes_inat'] = details.get('notes', 'Error')
            except Exception as e:
                print("UNABLE TO RETRIEVE INAT DATA FOR" + row['sample_code'] + f" ERROR: {e}")
    return samples







