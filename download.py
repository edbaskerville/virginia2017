import requests

def download_url(url, filename):
    done = False
    while not done:
        if os.path.exists(filename):
            size_so_far = os.path.getsize(filename)
        else:
            size_so_far = 0
        headers = {'Range': 'bytes={}-'.format(size_so_far)}
        try:
            response = requests.get(
                url, stream = True, headers = headers, timeout = 5.0
            )
            try:
                response_size = int(response.headers.get('Content-Length'))
            except:
                response_size = None
            if response_size is None:
                total_size = None
                size_so_far = 0
            else:
                total_size = size_so_far + response_size
            print(total_size)
            if response_size is None or response_size > 0:
                if size_so_far > 0:
                    print('Continuing incomplete download of {}...'.format(url))
                    mode = 'ab'
                else:
                    print('Starting download of {}...'.format(url))
                    mode = 'wb'
                
                with open(filename, mode) as file:
                    for block in response.iter_content(1024 * 128):
                        file.write(block)
                        size_so_far += len(block)
                        
                        size_so_far_mb = size_so_far / 1024.0 / 1024.0
                        if total_size is None:
                            sys.stdout.write('\r{:.2f} MB downloaded'.format(size_so_far_mb))
                        else:
                            total_size_mb = total_size / 1024.0 / 1024.0
                            percent_done = size_so_far_mb / total_size_mb * 100
                            sys.stdout.write(
                                '\r{:.2f} MB / {:.2f} MB downloaded ({:.2f}%)'.format(
                                    size_so_far_mb, total_size_mb, percent_done
                                )
                            )
                    print('')
                if total_size is not None:
                    assert size_so_far == total_size
                done = True
                
                if response.status_code != 200:
                    return response.status_code
            else:
                print('Download of {}) already finished.'.format(url))
                done = True
        except requests.exceptions.ConnectionError:
            print('')
            print('Exception occurred; retrying in 5 seconds...')
            time.sleep(5)
