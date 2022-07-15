def config_template(config):
    with open(config, 'r', encoding='utf-8') as f:
        config_dict = {}
        for line in f:
            if line.strip().startswith('#') or line.strip() == '':
                continue
            k, v = line.strip().split()
            config_dict[k] = str(v).strip()
    return config_dict

def config_template2(config):
    with open(config, 'r', encoding='utf-8') as f:
        config_dict = {}
        for line in f:
            if line.strip().startswith('#') or line.strip() == '':
                continue
            k, v = line.strip().split()
            config_dict[k] = str(v).strip()
    return config_dict

print(config_template('/Users/jmzhang/workspace/mygithub/bioinfor_utils/others/config.cfg'))
print(config_template2('/Users/jmzhang/workspace/mygithub/bioinfor_utils/others/config.cfg'))