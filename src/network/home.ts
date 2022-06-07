import {ins} from './request'

// 查询全部
export function uselist(){
  return ins().get('enum/?format=api')
}
// 增加操作，需要从实例中传入需要添加的数据
export function usepost(title:string,core:object,ligand:object){
  return ins().post('enum/?format=api',{
    title:title,
    core:core,
		ligand:ligand
  })
}
// 修改操作，需要从实例中传入需要修改的数据
export function useput(id:number,core:object,ligand:object){
  return ins().put('enum_mol/?format=api'+String(id),{
			core:core,
			ligand:ligand
    }
  )
}
// 删除操作，需要从实例中传入id
export function usedelete(id:number){
  return ins().delete('mol_design/enum_mol/?format=api'+String(id))
}
