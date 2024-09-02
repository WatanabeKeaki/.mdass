def ignore_comment(line):
    config_dict = {}
    line = line.strip()
    
    # '#' でコメントを切り捨て、両端の空白を削除
    line = line.split('#')[0].strip()
    if line:  # 空行でないか確認
        # '=' でキーと値を分割
        key_value = line.split('=')
        if len(key_value) == 2:  # 有効なキー=値のペア
            key = key_value[0].strip()
            value = key_value[1].strip()
            
            config_dict[key] = value
    
    return config_dict